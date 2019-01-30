module Numeric.LAPACK.Linear.LowerUpper (
   LowerUpper,
   Square,
   Transposition(..),
   Conjugation(..),
   Inversion(..),
   mapExtent,
   fromMatrix,
   toMatrix,
   solve,
   multiplyFullRight,

   determinant,

   extractP,
   multiplyP,

   extractL,
   wideExtractL,
   wideMultiplyL,
   wideSolveL,

   extractU,
   tallExtractU,
   tallMultiplyU,
   tallSolveU,

   caseTallWide,
   ) where

import qualified Numeric.LAPACK.Matrix.Multiply as MM
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Matrix.Basic as Basic
import qualified Numeric.LAPACK.Matrix.Private as Matrix
import qualified Numeric.LAPACK.Permutation.Private as Perm
import qualified Numeric.LAPACK.Split as Split
import Numeric.LAPACK.Matrix.Triangular.Basic (UnitLower, Upper)
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor, ColumnMajor), Triangle(Triangle))
import Numeric.LAPACK.Matrix.Private
         (Full, ZeroInt, zeroInt,
          Transposition(NonTransposed, Transposed),
          Conjugation(NonConjugated, Conjugated),
          Inversion(NonInverted, Inverted), flipInversion)
import Numeric.LAPACK.Linear.Private (solver, withInfo)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Format (Format(format))
import Numeric.LAPACK.Private
         (pointerSeq, peekCInt,
          copyBlock, copyTransposed, copyToColumnMajor)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked.Monadic as ArrayIO
import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.Marshal.Array (advancePtr)
import Foreign.C.Types (CInt)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Ptr (Ptr)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Monad (forM_)
import Control.Applicative ((<$>))


data LowerUpper vert horiz height width a =
   LowerUpper {
      _pivot :: Vector ZeroInt CInt,
      split_ ::
         Array
            (MatrixShape.Split MatrixShape.Triangle vert horiz height width) a
   } deriving (Show)

type Square sh = LowerUpper Extent.Small Extent.Small sh sh

instance
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
      Format (LowerUpper vert horiz height width a) where
   format fmt lu@(LowerUpper _ipiv m) = format fmt (extractP NonInverted lu, m)

mapExtent ::
   (Extent.C vertA, Extent.C horizA) =>
   (Extent.C vertB, Extent.C horizB) =>
   Extent.Map vertA horizA vertB horizB height width ->
   LowerUpper vertA horizA height width a ->
   LowerUpper vertB horizB height width a
mapExtent f (LowerUpper pivot split) =
   LowerUpper pivot $ Array.mapShape (MatrixShape.splitMapExtent f) split

{- |
@LowerUpper.fromMatrix a@
computes the LU decomposition of matrix @a@ with row pivotisation.

You can reconstruct @a@ from @lu@ depending on wether @a@ is tall or wide.

> LU.multiplyP False lu $ LU.extractL lu <#> LU.tallExtractU lu
> LU.multiplyP False lu $ LU.wideExtractL lu <#> LU.extractU lu
-}
fromMatrix ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Full vert horiz height width a ->
   LowerUpper vert horiz height width a
fromMatrix (Array (MatrixShape.Full order extent) a) =
   let (height,width) = Extent.dimensions extent
       m = Shape.size height
       n = Shape.size width
   in uncurry LowerUpper $
      Array.unsafeCreateWithSizeAndResult (zeroInt $ min m n) $ \_ ipivPtr ->
      ArrayIO.unsafeCreate
         (MatrixShape.Split MatrixShape.Triangle ColumnMajor extent) $ \luPtr ->

   evalContT $ do
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim m
      liftIO $ do
         copyToColumnMajor order m n aPtr luPtr
         withInfo "getrf" $ LapackGen.getrf mPtr nPtr luPtr ldaPtr ipivPtr

solve ::
   (Extent.C vert, Extent.C horiz, Eq height, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Square height a ->
   Full vert horiz height width a ->
   Full vert horiz height width a
solve
   (LowerUpper
      (Array _ ipiv)
      (Array (MatrixShape.Split MatrixShape.Triangle orderLU extentLU) lu)) =

   solver "LowerUpper.solve" (Extent.squareSize extentLU) $
         \n nPtr nrhsPtr xPtr ldxPtr -> do
      let lda = n
      transPtr <- Call.char 'N'
      aPtr <-
         case orderLU of
            RowMajor -> do
               aPtr <- ContT $ withForeignPtr lu
               atmpPtr <- Call.allocaArray (n*n)
               liftIO $ copyToColumnMajor orderLU n n aPtr atmpPtr
               return atmpPtr
            ColumnMajor -> ContT $ withForeignPtr lu
      ldaPtr <- Call.leadingDim lda
      ipivPtr <- ContT $ withForeignPtr ipiv
      liftIO $
         withInfo "getrs" $
            LapackGen.getrs transPtr
               nPtr nrhsPtr aPtr ldaPtr ipivPtr xPtr ldxPtr

{- |
Caution:
@LU.determinant . LU.fromMatrix@ will fail for singular matrices.
-}
determinant :: (Shape.C sh, Class.Floating a, Eq a) => Square sh a -> a
determinant (LowerUpper ipiv split) =
   let det = Split.determinantR split
   in if Split.oddPermutation $ Array.toList ipiv then -det else det


extractP ::
   (Extent.C vert, Extent.C horiz, Shape.C height) =>
   Inversion -> LowerUpper vert horiz height width a -> Perm.Permutation height
extractP inverted (LowerUpper ipiv (Array shape _)) =
   Perm.fromPivots (flipInversion inverted) (MatrixShape.splitHeight shape) ipiv

multiplyP ::
   (Extent.C vertA, Extent.C horizA, Extent.C vertB, Extent.C horizB,
    Eq height, Shape.C height, Shape.C widthA, Shape.C widthB,
    Class.Floating a) =>
   Inversion ->
   LowerUpper vertA horizA height widthA a ->
   Full vertB horizB height widthB a ->
   Full vertB horizB height widthB a
multiplyP inverted
      (LowerUpper (Array shapeIPiv ipiv)
         (Array (MatrixShape.Split _ _ extentLU) _lu))
      (Array shape@(MatrixShape.Full order extent) a) =
   Array.unsafeCreate shape $ \bPtr -> do

   Call.assert "multiplyP: heights mismatch"
      (Extent.height extentLU == Extent.height extent)

   let (height,width) = Extent.dimensions extent
   let m = Shape.size height
   let n = Shape.size width
   let k = Shape.size shapeIPiv

   evalContT $ do
      aPtr <- ContT $ withForeignPtr a
      ipivPtr <- ContT $ withForeignPtr ipiv
      liftIO $ copyBlock (n*m) aPtr bPtr
      case order of
         ColumnMajor -> do
            nPtr <- Call.cint n
            ldaPtr <- Call.leadingDim m
            k1Ptr <- Call.cint 1
            k2Ptr <- Call.cint k
            incxPtr <-
               Call.cint $
               case inverted of
                  Inverted -> 1
                  NonInverted -> -1
            liftIO $
               LapackGen.laswp nPtr bPtr ldaPtr k1Ptr k2Ptr ipivPtr incxPtr
         RowMajor ->
            liftIO $ swapColumns m bPtr $ take k $
            case inverted of
               Inverted -> zip [0..] $ pointerSeq 1 ipivPtr
               NonInverted ->
                  zip (iterate (subtract 1) (k-1)) $
                  pointerSeq (-1) (advancePtr ipivPtr (k-1))

{-# INLINE swapColumns #-}
swapColumns ::
   (Class.Floating a) =>
   Int -> Ptr a -> [(Int, Ptr CInt)] -> IO ()
swapColumns m xPtr ptrs = evalContT $ do
   mPtr <- Call.cint m
   incPtr <- Call.cint 1
   let columnPtr k = advancePtr xPtr (m*k)
   liftIO $ forM_ ptrs $ \(i,ipivPtr) -> do
      j <- subtract 1 <$> peekCInt ipivPtr
      BlasGen.swap mPtr (columnPtr i) incPtr (columnPtr j) incPtr



extractL ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   LowerUpper vert horiz height width a ->
   Full vert horiz height width a
extractL = Split.extractTriangle (Left Triangle) . split_

wideExtractL ::
   (Extent.C horiz, Shape.C height, Shape.C width, Class.Floating a) =>
   LowerUpper Extent.Small horiz height width a -> UnitLower height a
wideExtractL = Split.wideExtractL . split_

{- |
@wideMultiplyL transposed lu a@ multiplies the square part of @lu@
containing the lower triangular matrix with @a@.

> wideMultiplyL False lu a == wideExtractL lu <#> a
> wideMultiplyL True lu a == wideExtractL (Tri.transposeUp lu) <#> a
-}
wideMultiplyL ::
   (Extent.C horizA, Extent.C vert, Extent.C horiz, Shape.C height, Eq height,
    Shape.C widthA, Shape.C widthB, Class.Floating a) =>
   Transposition ->
   LowerUpper Extent.Small horizA height widthA a ->
   Full vert horiz height widthB a ->
   Full vert horiz height widthB a
wideMultiplyL transposed = Split.wideMultiplyL transposed . split_

wideSolveL ::
   (Extent.C horizA, Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Shape.C nrhs, Class.Floating a) =>
   Transposition -> Conjugation ->
   LowerUpper Extent.Small horizA height width a ->
   Full vert horiz height nrhs a -> Full vert horiz height nrhs a
wideSolveL transposed conjugated =
   Split.wideSolveL transposed conjugated . split_


extractU ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   LowerUpper vert horiz height width a ->
   Full vert horiz height width a
extractU = Split.extractTriangle (Right Triangle) . split_

tallExtractU ::
   (Extent.C vert, Shape.C height, Shape.C width, Class.Floating a) =>
   LowerUpper vert Extent.Small height width a -> Upper width a
tallExtractU = Split.tallExtractR . split_

{- |
@tallMultiplyU transposed lu a@ multiplies the square part of @lu@
containing the upper triangular matrix with @a@.

> tallMultiplyU False lu a == tallExtractU lu <#> a
> tallMultiplyU True lu a == tallExtractU (Tri.transposeDown lu) <#> a
-}
tallMultiplyU ::
   (Extent.C vertA, Extent.C vert, Extent.C horiz, Shape.C height, Eq height,
    Shape.C heightA, Shape.C widthB, Class.Floating a) =>
   Transposition ->
   LowerUpper vertA Extent.Small heightA height a ->
   Full vert horiz height widthB a ->
   Full vert horiz height widthB a
tallMultiplyU transposed = Split.tallMultiplyR transposed . split_

tallSolveU ::
   (Extent.C vertA, Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Eq width, Shape.C nrhs, Class.Floating a) =>
   Transposition -> Conjugation ->
   LowerUpper vertA Extent.Small height width a ->
   Full vert horiz width nrhs a -> Full vert horiz width nrhs a
tallSolveU transposed conjugated =
   Split.tallSolveR transposed conjugated . split_



toMatrix ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Eq width, Class.Floating a) =>
   LowerUpper vert horiz height width a ->
   Full vert horiz height width a
toMatrix =
   getToMatrix $
   Extent.switchTagPair
      (ToMatrix wideToMatrix)
      (ToMatrix wideToMatrix)
      (ToMatrix tallToMatrix)
      (ToMatrix $
         either
            (Matrix.fromFull . tallToMatrix)
            (Matrix.fromFull . wideToMatrix) .
         caseTallWide)

newtype ToMatrix height width a vert horiz =
   ToMatrix {
      getToMatrix ::
         LowerUpper vert horiz height width a ->
         Full vert horiz height width a
   }

tallToMatrix ::
   (Extent.C vert, Shape.C height, Shape.C width, Eq height, Eq width,
    Class.Floating a) =>
   LowerUpper vert Extent.Small height width a ->
   Full vert Extent.Small height width a
tallToMatrix a =
   multiplyP NonInverted a $ Basic.transpose $
   tallMultiplyU Transposed a $ Basic.transpose $ extractL a

wideToMatrix ::
   (Extent.C horiz, Shape.C height, Shape.C width, Eq height, Eq width,
    Class.Floating a) =>
   LowerUpper Extent.Small horiz height width a ->
   Full Extent.Small horiz height width a
wideToMatrix a =
   multiplyP NonInverted a $ wideMultiplyL NonTransposed a $ extractU a


multiplyFullRight ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Shape.C fuse, Eq fuse,
    Class.Floating a) =>
   LowerUpper vert horiz height fuse a ->
   Full vert horiz fuse width a ->
   Full vert horiz height width a
multiplyFullRight =
   getMultiplyFullRight $
   Extent.switchTagPair
      (MultiplyFullRight wideMultiplyFullRight)
      (MultiplyFullRight wideMultiplyFullRight)
      (MultiplyFullRight tallMultiplyFullRight)
      (MultiplyFullRight $
         either tallMultiplyFullRight wideMultiplyFullRight . caseTallWide)

newtype MultiplyFullRight height fuse width a vert horiz =
   MultiplyFullRight {
      getMultiplyFullRight ::
         LowerUpper vert horiz height fuse a ->
         Full vert horiz fuse width a ->
         Full vert horiz height width a
   }

tallMultiplyFullRight ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Shape.C fuse, Eq height, Eq fuse,
    Class.Floating a) =>
   LowerUpper vert Extent.Small height fuse a ->
   Full vert horiz fuse width a ->
   Full vert horiz height width a
tallMultiplyFullRight a =
   multiplyP NonInverted a .
   MM.multiply (Matrix.generalizeTall (extractL a)) .
   tallMultiplyU NonTransposed a

wideMultiplyFullRight ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Shape.C fuse, Eq height, Eq fuse,
    Class.Floating a) =>
   LowerUpper Extent.Small horiz height fuse a ->
   Full vert horiz fuse width a ->
   Full vert horiz height width a
wideMultiplyFullRight a =
   multiplyP NonInverted a . wideMultiplyL NonTransposed a .
   MM.multiply (Matrix.generalizeWide (extractU a))


type Tall = LowerUpper Extent.Big Extent.Small
type Wide = LowerUpper Extent.Small Extent.Big

caseTallWide ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
   LowerUpper vert horiz height width a ->
   Either (Tall height width a) (Wide height width a)
caseTallWide (LowerUpper ipiv (Array shape a)) =
   either
      (Left . LowerUpper ipiv . flip Array a)
      (Right . LowerUpper ipiv . flip Array a) $
   MatrixShape.caseTallWideSplit shape


_toRowMajor ::
   (Extent.C vert, Extent.C horiz, Eq height, Shape.C height, Shape.C width,
    Class.Floating a) =>
   LowerUpper vert horiz height width a ->
   LowerUpper vert horiz height width a
_toRowMajor
   (LowerUpper ipiv
      arr@(Array (MatrixShape.Split MatrixShape.Triangle order extent) a)) =
   LowerUpper ipiv $
   case order of
      RowMajor -> arr
      ColumnMajor ->
         Array.unsafeCreate
            (MatrixShape.Split MatrixShape.Triangle RowMajor extent) $ \bPtr ->
         withForeignPtr a $ \aPtr -> do
            let (height, width) = Extent.dimensions extent
            let n = Shape.size width
            let m = Shape.size height
            copyTransposed n m aPtr n bPtr
