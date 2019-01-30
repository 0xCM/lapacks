{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Numeric.LAPACK.Matrix (
   Full,
   General, Tall, Wide,
   ZeroInt, zeroInt,
   transpose, adjoint,
   Matrix.height, Matrix.width,
   caseTallWide,
   fromScalar, toScalar,
   fromList,
   mapExtent, fromFull,
   generalizeTall, generalizeWide,
   identity,
   diagonal,
   fromRowsNonEmpty,    fromRowArray,    fromRows,
   fromColumnsNonEmpty, fromColumnArray, fromColumns,
   Basic.singleRow,   Basic.singleColumn,
   Basic.flattenRow,  Basic.flattenColumn,
   toRows, toColumns,
   toRowArray, toColumnArray,
   takeRow, takeColumn,
   takeRows, takeColumns, takeEqually,
   dropRows, dropColumns, dropEqually,
   takeTopRows, takeBottomRows,
   takeLeftColumns, takeRightColumns,
   reverseRows, reverseColumns,
   fromRowMajor, toRowMajor, flatten,
   forceOrder, adaptOrder,
   (|||),
   (===),

   tensorProduct,
   outer,
   sumRank1,

   RealOf,
   add, sub,
   rowSums, columnSums,
   scaleRows, scaleColumns,
   scaleRowsComplex, scaleColumnsComplex,
   scaleRowsReal, scaleColumnsReal,
   multiply,
   multiplyVector,

   Multiply, (<#>),
   MultiplyLeft, (<#),
   MultiplyRight, (#>),

   Solve, solve, solveVector,
   Inverse, inverse,
   ) where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Square.Basic as Square
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Matrix.Basic as Basic
import qualified Numeric.LAPACK.Matrix.Private as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Matrix.Shape.Private (Order(RowMajor, ColumnMajor))
import Numeric.LAPACK.Matrix.Multiply
         (Multiply((<#>)), MultiplyLeft((<#)), MultiplyRight((#>)),
          multiplyVector, multiply, multiplyVectorUnchecked)
import Numeric.LAPACK.Matrix.Divide
         (Solve(solve), solveVector, Inverse(inverse))
import Numeric.LAPACK.Matrix.Basic
         (transpose, forceOrder, forceRowMajor, scaleRows, scaleColumns)
import Numeric.LAPACK.Matrix.Private
         (Full, Tall, Wide, General, argGeneral, ZeroInt, zeroInt,
          mapExtent, fromFull, generalizeTall, generalizeWide)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, zero, one)
import Numeric.LAPACK.Private
         (pointerSeq, fill, copyTransposed, copySubMatrix, copyBlock)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Boxed as BoxedArray
import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))
import Data.Array.Comfort.Shape ((:+:)((:+:)))

import Foreign.Marshal.Array (copyArray, advancePtr, pokeArray)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Ptr (Ptr, castPtr)
import Foreign.Storable (Storable, poke, peek)

import System.IO.Unsafe (unsafePerformIO)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)

import qualified Data.NonEmpty as NonEmpty
import Data.Complex (Complex)
import Data.Foldable (forM_)
import Data.Bool.HT (if')


{- |
conjugate transpose

Problem: @adjoint a <#> a@ is always square,
but how to convince the type checker to choose the Square type?
Anser: Use @Hermitian.toSquare $ Hermitian.covariance a@ instead.
-}
adjoint ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Full vert horiz height width a -> Full horiz vert width height a
adjoint = transpose . Vector.conjugate


{- |
Square matrices will be classified as 'Tall'.
-}
caseTallWide ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
   Full vert horiz height width a ->
   Either (Tall height width a) (Wide height width a)
caseTallWide (Array shape a) =
   either (Left . flip Array a) (Right . flip Array a) $
   MatrixShape.caseTallWide shape


fromScalar :: (Storable a) => a -> General () () a
fromScalar = Square.toGeneral . Square.fromScalar

toScalar :: (Storable a) => General () () a -> a
toScalar = argGeneral $ \_ () () a ->
   unsafePerformIO $ withForeignPtr a peek

fromList ::
   (Shape.C height, Shape.C width, Storable a) =>
   height -> width -> [a] -> General height width a
fromList height width =
   Array.fromList (MatrixShape.general RowMajor height width)


identity ::
   (Shape.C sh, Class.Floating a) =>
   sh -> General sh sh a
identity = Square.toGeneral . Square.identity

diagonal ::
   (Shape.C sh, Class.Floating a) =>
   Vector sh a -> General sh sh a
diagonal = Square.toGeneral . Square.diagonal


fromRowsNonEmpty ::
   (Shape.C width, Eq width, Storable a) =>
   NonEmpty.T [] (Vector width a) -> General ZeroInt width a
fromRowsNonEmpty (NonEmpty.Cons row rows) =
   fromRows (Array.shape row) (row:rows)

fromRowArray ::
   (Shape.C height, Shape.C width, Eq width, Storable a) =>
   width -> BoxedArray.Array height (Vector width a) -> General height width a
fromRowArray width rows =
   Array.reshape (MatrixShape.general RowMajor (BoxedArray.shape rows) width) $
   fromRows width $ BoxedArray.toList rows

fromRows ::
   (Shape.C width, Eq width, Storable a) =>
   width -> [Vector width a] -> General ZeroInt width a
fromRows width rows =
   Array.unsafeCreate
      (MatrixShape.general RowMajor (zeroInt $ length rows) width)
      (gather width rows)

fromColumnsNonEmpty ::
   (Shape.C height, Eq height, Storable a) =>
   NonEmpty.T [] (Vector height a) -> General height ZeroInt a
fromColumnsNonEmpty (NonEmpty.Cons column columns) =
   fromColumns (Array.shape column) (column:columns)

fromColumnArray ::
   (Shape.C height, Eq height, Shape.C width, Storable a) =>
   height -> BoxedArray.Array width (Vector height a) -> General height width a
fromColumnArray height columns =
   Array.reshape
      (MatrixShape.general ColumnMajor height (BoxedArray.shape columns)) $
   fromColumns height $ BoxedArray.toList columns

fromColumns ::
   (Shape.C height, Eq height, Storable a) =>
   height -> [Vector height a] -> General height ZeroInt a
fromColumns height columns =
   Array.unsafeCreate
      (MatrixShape.general ColumnMajor height (zeroInt $ length columns))
      (gather height columns)

gather ::
   (Shape.C width, Eq width, Storable a) =>
   width -> [Array width a] -> Ptr a -> IO ()
gather width rows dstPtr =
   let widthSize = Shape.size width
   in forM_ (zip (pointerSeq widthSize dstPtr) rows) $
         \(dstRowPtr, Array.Array rowWidth srcFPtr) ->
         withForeignPtr srcFPtr $ \srcPtr -> do
            Call.assert
               "Matrix.fromRows/fromColumnsNonEmpty: non-matching vector size"
               (width == rowWidth)
            copyArray dstRowPtr srcPtr widthSize


toRows ::
   (Extent.C vert, Extent.C horiz,
    Shape.Indexed height, Shape.C width, Class.Floating a) =>
   Full vert horiz height width a -> [Vector width a]
toRows a = map (takeRow a) $ Shape.indices $ Matrix.height a

toColumns ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.Indexed width, Class.Floating a) =>
   Full vert horiz height width a -> [Vector height a]
toColumns a = map (takeColumn a) $ Shape.indices $ Matrix.width a

toRowArray ::
   (Extent.C vert, Extent.C horiz,
    Shape.Indexed height, Shape.C width, Class.Floating a) =>
   Full vert horiz height width a -> BoxedArray.Array height (Vector width a)
toRowArray a =
   let height = Matrix.height a
   in BoxedArray.fromList height $ map (takeRow a) $ Shape.indices height

toColumnArray ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.Indexed width, Class.Floating a) =>
   Full vert horiz height width a -> BoxedArray.Array width (Vector height a)
toColumnArray a =
   let width = Matrix.width a
   in BoxedArray.fromList width $ map (takeColumn a) $ Shape.indices width


takeRow ::
   (Extent.C vert, Extent.C horiz,
    Shape.Indexed height, Shape.C width, Shape.Index height ~ ix,
    Class.Floating a) =>
   Full vert horiz height width a -> ix -> Vector width a
takeRow (Array (MatrixShape.Full order extent) x) ix =
   let (height,width) = Extent.dimensions extent
   in case order of
         RowMajor -> pickConsecutive height width x ix
         ColumnMajor -> pickScattered width height x ix

takeColumn ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.Indexed width, Shape.Index width ~ ix,
    Class.Floating a) =>
   Full vert horiz height width a -> ix -> Vector height a
takeColumn (Array (MatrixShape.Full order extent) x) ix =
   let (height,width) = Extent.dimensions extent
   in case order of
         RowMajor -> pickScattered height width x ix
         ColumnMajor -> pickConsecutive width height x ix

pickConsecutive ::
   (Shape.Indexed height, Shape.C width, Shape.Index height ~ ix,
    Class.Floating a) =>
   height -> width -> ForeignPtr a -> ix -> Vector width a
pickConsecutive height width x ix =
   Array.unsafeCreateWithSize width $ \n yPtr -> evalContT $ do
      let offset = Shape.offset height ix
      nPtr <- Call.cint n
      xPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 1
      liftIO $
         BlasGen.copy nPtr (advancePtr xPtr (n*offset)) incxPtr yPtr incyPtr

pickScattered ::
   (Shape.C height, Shape.Indexed width, Shape.Index width ~ ix,
    Class.Floating a) =>
   height -> width -> ForeignPtr a -> ix -> Vector height a
pickScattered height width x ix =
   Array.unsafeCreateWithSize height $ \n yPtr -> evalContT $ do
      let offset = Shape.offset width ix
      nPtr <- Call.cint n
      xPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint $ Shape.size width
      incyPtr <- Call.cint 1
      liftIO $
         BlasGen.copy nPtr (advancePtr xPtr offset) incxPtr yPtr incyPtr


takeTopRows ::
   (Extent.C vert, Shape.C height0, Shape.C height1, Shape.C width,
    Class.Floating a) =>
   Full vert Extent.Big (height0:+:height1) width a ->
   Full vert Extent.Big height0 width a
takeTopRows (Array (MatrixShape.Full order extentA) a) =
   let (heightA@(heightB:+:_), width) = Extent.dimensions extentA
       extentB = Extent.reduceWideHeight heightB extentA
       ma = Shape.size heightA
       mb = Shape.size heightB
       n = Shape.size width
   in Array.unsafeCreateWithSize (MatrixShape.Full order extentB) $
            \blockSize bPtr ->
      withForeignPtr a $ \aPtr ->
      case order of
         RowMajor -> copyBlock blockSize aPtr bPtr
         ColumnMajor -> copySubMatrix mb n ma aPtr mb bPtr

takeBottomRows ::
   (Extent.C vert, Shape.C height0, Shape.C height1, Shape.C width,
    Class.Floating a) =>
   Full vert Extent.Big (height0:+:height1) width a ->
   Full vert Extent.Big height1 width a
takeBottomRows (Array (MatrixShape.Full order extentA) a) =
   let (heightA@(height0:+:heightB), width) = Extent.dimensions extentA
       extentB = Extent.reduceWideHeight heightB extentA
       k = Shape.size height0
       ma = Shape.size heightA
       mb = Shape.size heightB
       n = Shape.size width
   in Array.unsafeCreateWithSize (MatrixShape.Full order extentB) $
            \blockSize bPtr ->
      withForeignPtr a $ \aPtr ->
      case order of
         RowMajor -> copyBlock blockSize (advancePtr aPtr (k*n)) bPtr
         ColumnMajor -> copySubMatrix mb n ma (advancePtr aPtr k) mb bPtr

takeLeftColumns ::
   (Extent.C vert, Shape.C height, Shape.C width0, Shape.C width1,
    Class.Floating a) =>
   Full Extent.Big vert height (width0:+:width1) a ->
   Full Extent.Big vert height width0 a
takeLeftColumns = transpose . takeTopRows . transpose

takeRightColumns ::
   (Extent.C vert, Shape.C height, Shape.C width0, Shape.C width1,
    Class.Floating a) =>
   Full Extent.Big vert height (width0:+:width1) a ->
   Full Extent.Big vert height width1 a
takeRightColumns = transpose . takeBottomRows . transpose


splitRows ::
   (Extent.C vert, Shape.C width, Class.Floating a) =>
   Int ->
   Full vert Extent.Big ZeroInt width a ->
   Full vert Extent.Big (ZeroInt:+:ZeroInt) width a
splitRows k =
   Array.mapShape
      (\(MatrixShape.Full order extentA) ->
         let (Shape.ZeroBased heightA) = Extent.height extentA
             heightB = min k heightA
         in if' (k<0) (error "split: negative number") $
            MatrixShape.Full order $
            Extent.reduceWideHeight
               (Shape.ZeroBased heightB :+: Shape.ZeroBased (heightA-heightB))
               extentA)

takeRows, dropRows ::
   (Extent.C vert, Shape.C width, Class.Floating a) =>
   Int ->
   Full vert Extent.Big ZeroInt width a ->
   Full vert Extent.Big ZeroInt width a
takeRows k = takeTopRows . splitRows k
dropRows k = takeBottomRows . splitRows k

takeColumns, dropColumns ::
   (Extent.C horiz, Shape.C height, Class.Floating a) =>
   Int ->
   Full Extent.Big horiz height ZeroInt a ->
   Full Extent.Big horiz height ZeroInt a
takeColumns k = transpose . takeRows k . transpose
dropColumns k = transpose . dropRows k . transpose


{- |
Take a left-top aligned square or as much as possible of it.
The advantange of this function is that it maintains the matrix size relation,
e.g. Square remains Square, Tall remains Tall.
-}
takeEqually ::
   (Extent.C vert, Extent.C horiz, Class.Floating a) =>
   Int ->
   Full vert horiz ZeroInt ZeroInt a ->
   Full vert horiz ZeroInt ZeroInt a
takeEqually k (Array (MatrixShape.Full order extentA) a) =
   let (Shape.ZeroBased heightA, Shape.ZeroBased widthA) =
         Extent.dimensions extentA
       heightB = min k heightA
       widthB  = min k widthA
       extentB =
         Extent.reduceConsistent
            (Shape.ZeroBased heightB) (Shape.ZeroBased widthB) extentA
   in if' (k<0) (error "take: negative number") $
      Array.unsafeCreate (MatrixShape.Full order extentB) $ \bPtr ->
      withForeignPtr a $ \aPtr ->
      case order of
         RowMajor -> copySubMatrix widthB heightB widthA aPtr widthB bPtr
         ColumnMajor -> copySubMatrix heightB widthB heightA aPtr heightB bPtr

{- |
Drop the same number of top-most rows and left-most columns.
The advantange of this function is that it maintains the matrix size relation,
e.g. Square remains Square, Tall remains Tall.
-}
dropEqually ::
   (Extent.C vert, Extent.C horiz, Class.Floating a) =>
   Int ->
   Full vert horiz ZeroInt ZeroInt a ->
   Full vert horiz ZeroInt ZeroInt a
dropEqually k (Array (MatrixShape.Full order extentA) a) =
   let (Shape.ZeroBased heightA, Shape.ZeroBased widthA) =
         Extent.dimensions extentA
       heightB = heightA - top; top  = min k heightA
       widthB  = widthA - left; left = min k widthA
       extentB =
         Extent.reduceConsistent
            (Shape.ZeroBased heightB) (Shape.ZeroBased widthB) extentA
   in if' (k<0) (error "drop: negative number") $
      Array.unsafeCreate (MatrixShape.Full order extentB) $ \bPtr ->
      withForeignPtr a $ \aPtr ->
      case order of
         RowMajor ->
            copySubMatrix widthB heightB
               widthA (advancePtr aPtr (top*widthA+left)) widthB bPtr
         ColumnMajor ->
            copySubMatrix heightB widthB
               heightA (advancePtr aPtr (left*heightA+top)) heightB bPtr


-- alternative: laswp
reverseRows ::
   (Extent.C vert, Extent.C horiz, Shape.C width, Class.Floating a) =>
   Full vert horiz ZeroInt width a -> Full vert horiz ZeroInt width a
reverseRows (Array shape@(MatrixShape.Full order extent) a) =
   Array.unsafeCreateWithSize shape $ \blockSize bPtr -> evalContT $ do
      let (height,width) = Extent.dimensions extent
      let n = Shape.size height
      let m = Shape.size width
      fwdPtr <- Call.bool True
      nPtr <- Call.cint n
      mPtr <- Call.cint m
      kPtr <- Call.allocaArray n
      aPtr <- ContT $ withForeignPtr a
      liftIO $ do
         copyBlock blockSize aPtr bPtr
         pokeArray kPtr $ take n $ iterate (subtract 1) $ fromIntegral n
         case order of
            RowMajor -> LapackGen.lapmt fwdPtr mPtr nPtr bPtr mPtr kPtr
            ColumnMajor -> LapackGen.lapmr fwdPtr nPtr mPtr bPtr nPtr kPtr

reverseColumns ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Class.Floating a) =>
   Full vert horiz height ZeroInt a -> Full vert horiz height ZeroInt a
reverseColumns = transpose . reverseRows . transpose


fromRowMajor ::
   (Shape.C height, Shape.C width, Class.Floating a) =>
   Array (height,width) a -> General height width a
fromRowMajor = Array.mapShape (uncurry $ MatrixShape.general RowMajor)

toRowMajor ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Class.Floating a) =>
   Full vert horiz height width a -> Array (height,width) a
toRowMajor =
   Array.mapShape
      (\shape -> (MatrixShape.fullHeight shape, MatrixShape.fullWidth shape)) .
   forceRowMajor

{- |
@adaptOrder x y@ contains the data of @y@ with the layout of @x@.
-}
adaptOrder ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Full vert horiz height width a ->
   Full vert horiz height width a ->
   Full vert horiz height width a
adaptOrder x = forceOrder (MatrixShape.fullOrder $ Array.shape x)

flatten ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Class.Floating a) =>
   Full vert horiz height width a -> Vector ZeroInt a
flatten = Array.mapShape (zeroInt . Shape.size) . toRowMajor


infixl 3 |||
infixl 2 ===

(|||) ::
   (Extent.C vert, Shape.C height, Eq height, Shape.C widtha, Shape.C widthb,
    Class.Floating a) =>
   Full vert Extent.Big height widtha a ->
   Full vert Extent.Big height widthb a ->
   Full vert Extent.Big height (widtha:+:widthb) a
(|||)
      (Array (MatrixShape.Full orderA extentA) a)
      (Array (MatrixShape.Full orderB extentB) b) =
   let (heightA,widthA) = Extent.dimensions extentA
       (heightB,widthB) = Extent.dimensions extentB
       extent = Extent.widen (widthA:+:widthB) extentA
       shape order = MatrixShape.Full order extent
   in
    if heightA /= heightB
      then error "(|||): mismatching heights"
      else
         case (orderA,orderB) of
            (RowMajor,RowMajor) ->
               Array.unsafeCreate (shape RowMajor) $
               \cPtr -> evalContT $ do
                  let n = Shape.size heightA
                  let ma = Shape.size widthA
                  let mb = Shape.size widthB
                  let m = ma+mb
                  maPtr <- Call.cint ma
                  mbPtr <- Call.cint mb
                  aPtr <- ContT $ withForeignPtr a
                  bPtr <- ContT $ withForeignPtr b
                  incxPtr <- Call.cint 1
                  incyPtr <- Call.cint 1
                  liftIO $
                     sequence_ $ take n $
                     zipWith3
                        (\akPtr bkPtr ckPtr -> do
                           BlasGen.copy maPtr akPtr incxPtr ckPtr incyPtr
                           BlasGen.copy mbPtr bkPtr incxPtr
                              (ckPtr `advancePtr` ma) incyPtr)
                        (pointerSeq ma aPtr)
                        (pointerSeq mb bPtr)
                        (pointerSeq m cPtr)
            (RowMajor,ColumnMajor) ->
               Array.unsafeCreate (shape ColumnMajor) $
               \cPtr -> evalContT $ do
                  let n = Shape.size heightA
                  let ma = Shape.size widthA
                  let mb = Shape.size widthB
                  aPtr <- ContT $ withForeignPtr a
                  bPtr <- ContT $ withForeignPtr b
                  liftIO $ do
                     copyTransposed n ma aPtr n cPtr
                     copyBlock (n*mb) bPtr (advancePtr cPtr (n*ma))
            (ColumnMajor,RowMajor) ->
               Array.unsafeCreate (shape ColumnMajor) $
               \cPtr -> evalContT $ do
                  let n = Shape.size heightA
                  let ma = Shape.size widthA
                  let mb = Shape.size widthB
                  let volA = n*ma
                  aPtr <- ContT $ withForeignPtr a
                  bPtr <- ContT $ withForeignPtr b
                  liftIO $ do
                     copyBlock volA aPtr cPtr
                     copyTransposed n mb bPtr n (advancePtr cPtr volA)
            (ColumnMajor,ColumnMajor) ->
               Array.unsafeCreate (shape ColumnMajor) $
               \cPtr -> evalContT $ do
                  let n = Shape.size heightA
                  let na = n * Shape.size widthA
                  let nb = n * Shape.size widthB
                  naPtr <- Call.cint na
                  nbPtr <- Call.cint nb
                  aPtr <- ContT $ withForeignPtr a
                  bPtr <- ContT $ withForeignPtr b
                  incxPtr <- Call.cint 1
                  incyPtr <- Call.cint 1
                  liftIO $ do
                     BlasGen.copy naPtr aPtr incxPtr cPtr incyPtr
                     BlasGen.copy nbPtr bPtr incxPtr
                        (cPtr `advancePtr` na) incyPtr

(===) ::
   (Extent.C horiz, Shape.C width, Eq width, Shape.C heighta, Shape.C heightb,
    Class.Floating a) =>
   Full Extent.Big horiz heighta width a ->
   Full Extent.Big horiz heightb width a ->
   Full Extent.Big horiz (heighta:+:heightb) width a
(===) a b = transpose (transpose a ||| transpose b)


add, sub ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Eq height, Eq width,
    Class.Floating a) =>
   Full vert horiz height width a ->
   Full vert horiz height width a ->
   Full vert horiz height width a
add x y = Vector.add (adaptOrder y x) y
sub x y = Vector.sub (adaptOrder y x) y


rowSums ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Class.Floating a) =>
   Full vert horiz height width a -> Vector height a
rowSums m =
   let width = MatrixShape.fullWidth $ Array.shape m
   in  multiplyVectorUnchecked m (Vector.constant width one)

columnSums ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Class.Floating a) =>
   Full vert horiz height width a -> Vector width a
columnSums m =
   let height = MatrixShape.fullHeight $ Array.shape m
   in  multiplyVectorUnchecked (transpose m) (Vector.constant height one)


scaleRowsComplex ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Class.Real a) =>
   Vector height a ->
   Full vert horiz height width (Complex a) ->
   Full vert horiz height width (Complex a)
scaleRowsComplex
   (Array heightX x) (Array shape@(MatrixShape.Full order extent) a) =
      Array.unsafeCreate shape $ \bComplexPtr -> do
   let (height,width) = Extent.dimensions extent
   Call.assert "scaleRowsComplex: sizes mismatch" (heightX == height)
   let bPtr = castPtr bComplexPtr
   case order of
      RowMajor -> evalContT $ do
         let m = Shape.size height
         let n = Shape.size width * 2
         alphaPtr <- Call.alloca
         nPtr <- Call.cint n
         xPtr <- ContT $ withForeignPtr x
         aPtr <- fmap castPtr $ ContT $ withForeignPtr a
         incaPtr <- Call.cint 1
         incbPtr <- Call.cint 1
         liftIO $ sequence_ $ take m $
            zipWith3
               (\xkPtr akPtr bkPtr -> do
                  poke alphaPtr =<< peek xkPtr
                  BlasGen.copy nPtr akPtr incaPtr bkPtr incbPtr
                  BlasGen.scal nPtr alphaPtr bkPtr incbPtr)
               (pointerSeq 1 xPtr)
               (pointerSeq n aPtr)
               (pointerSeq n bPtr)
      ColumnMajor -> evalContT $ do
         let m = Shape.size width
         let nr = Shape.size height
         let n = 2*nr
         transPtr <- Call.char 'N'
         nrPtr <- Call.cint nr
         nPtr <- Call.cint n
         klPtr <- Call.cint 0
         kuPtr <- Call.cint 0
         alphaPtr <- Call.number one
         xrPtr <- ContT $ withForeignPtr x
         xPtr <- Call.allocaArray n
         incxrPtr <- Call.cint 1
         incxPtr <- Call.cint 2
         ldxPtr <- Call.leadingDim 1
         aPtr <- fmap castPtr $ ContT $ withForeignPtr a
         incaPtr <- Call.cint 1
         betaPtr <- Call.number zero
         incbPtr <- Call.cint 1
         liftIO $ do
            BlasGen.copy nrPtr xrPtr incxrPtr xPtr incxPtr
            BlasGen.copy nrPtr xrPtr incxrPtr (advancePtr xPtr 1) incxPtr
            sequence_ $ take m $
               zipWith
                  (\akPtr bkPtr ->
                     Private.gbmv transPtr
                        nPtr nPtr klPtr kuPtr alphaPtr xPtr ldxPtr
                        akPtr incaPtr betaPtr bkPtr incbPtr)
                  (pointerSeq n aPtr)
                  (pointerSeq n bPtr)

scaleColumnsComplex ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Eq width, Class.Real a) =>
   Vector width a ->
   Full vert horiz height width (Complex a) ->
   Full vert horiz height width (Complex a)
scaleColumnsComplex x = transpose . scaleRowsComplex x . transpose


scaleRowsReal ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Eq height, Shape.C width,
    Class.Floating a) =>
   Vector height (RealOf a) ->
   Full vert horiz height width a ->
   Full vert horiz height width a
scaleRowsReal =
   getScaleRowsReal $
   Class.switchFloating
      (ScaleRowsReal scaleRows)
      (ScaleRowsReal scaleRows)
      (ScaleRowsReal scaleRowsComplex)
      (ScaleRowsReal scaleRowsComplex)

newtype ScaleRowsReal f g a =
   ScaleRowsReal {getScaleRowsReal :: f (RealOf a) -> g a -> g a}

scaleColumnsReal ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Eq width, Class.Floating a) =>
   Vector width (RealOf a) ->
   Full vert horiz height width a ->
   Full vert horiz height width a
scaleColumnsReal x = transpose . scaleRowsReal x . transpose


{- |
> tensorProduct order x y = singleColumn order x <#> singleRow order y
-}
tensorProduct ::
   (Shape.C height, Eq height, Shape.C width, Eq width, Class.Floating a) =>
   Order -> Vector height a -> Vector width a -> General height width a
tensorProduct order x y =
   case order of
      ColumnMajor -> tensorProd 'T' order x y
      RowMajor -> transpose $ tensorProd 'T' order y x

{- |
> outer order x y = tensorProduct order x (Vector.conjugate y)
-}
outer ::
   (Shape.C height, Eq height, Shape.C width, Eq width, Class.Floating a) =>
   Order -> Vector height a -> Vector width a -> General height width a
outer order x y =
   case order of
      ColumnMajor -> tensorProd 'C' ColumnMajor x y
      RowMajor -> transpose $ tensorProd 'C' RowMajor y x

{-# INLINE tensorProd #-}
tensorProd ::
   (Shape.C height, Eq height, Shape.C width, Eq width, Class.Floating a) =>
   Char -> Order ->
   Vector height a -> Vector width a -> General height width a
tensorProd trans order (Array shX x) (Array shY y) =
   Array.unsafeCreate (MatrixShape.general MatrixShape.ColumnMajor shX shY) $
      \cPtr -> do
   let m = Shape.size shX
   let n = Shape.size shY
   let ((transa,transb),(lda,ldb)) =
         case order of
            ColumnMajor -> (('N',trans),(m,n))
            RowMajor -> ((trans,'N'),(1,1))
   evalContT $ do
      transaPtr <- Call.char transa
      transbPtr <- Call.char transb
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      kPtr <- Call.cint 1
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr x
      ldaPtr <- Call.leadingDim lda
      bPtr <- ContT $ withForeignPtr y
      ldbPtr <- Call.leadingDim ldb
      betaPtr <- Call.number zero
      ldcPtr <- Call.leadingDim m
      liftIO $
         BlasGen.gemm
            transaPtr transbPtr mPtr nPtr kPtr alphaPtr
            aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr


sumRank1 ::
   (Shape.C height, Eq height, Shape.C width, Eq width, Class.Floating a) =>
   (height,width) ->
   [(a, (Vector height a, Vector width a))] -> General height width a
sumRank1 (height,width) xys =
   Array.unsafeCreateWithSize (MatrixShape.general ColumnMajor height width) $
      \size aPtr ->
   evalContT $ do
      let m = Shape.size height
      let n = Shape.size width
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.alloca
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 1
      ldaPtr <- Call.leadingDim m
      liftIO $ do
         fill zero size aPtr
         forM_ xys $ \(alpha, (Array shX x, Array shY y)) ->
            withForeignPtr x $ \xPtr ->
            withForeignPtr y $ \yPtr -> do
               Call.assert "Matrix.sumRank1: non-matching height" (height==shX)
               Call.assert "Matrix.sumRank1: non-matching width" (width==shY)
               poke alphaPtr alpha
               BlasGen.gerc mPtr nPtr
                  alphaPtr xPtr incxPtr yPtr incyPtr aPtr ldaPtr
