{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Singular (
   values,
   valuesTall,
   valuesWide,
   decompose,
   decomposeTall,
   decomposeWide,
   determinantAbsolute,
   leastSquaresMinimumNormRCond,
   pseudoInverseRCond,
   RealOf,
   ) where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Square.Basic as Square
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Matrix.Hermitian.Private
         (TakeDiagonal(..), Determinant(..))
import Numeric.LAPACK.Matrix.Extent.Private (Extent)
import Numeric.LAPACK.Matrix.Square.Basic (Square)
import Numeric.LAPACK.Matrix.Shape.Private (Order(ColumnMajor), swapOnRowMajor)
import Numeric.LAPACK.Matrix (scaleRowsReal)
import Numeric.LAPACK.Matrix.Private (Full, General, ZeroInt, zeroInt)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, zero)
import Numeric.LAPACK.Private
         (withAutoWorkspace, peekCInt, createHigherArray,
          copyToTemp, copyToColumnMajor, copyToSubColumnMajor)

import qualified Numeric.LAPACK.FFI.Complex as LapackComplex
import qualified Numeric.LAPACK.FFI.Real as LapackReal
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked.Monadic as ArrayIO
import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import System.IO.Unsafe (unsafePerformIO)

import qualified Foreign.Marshal.Array.Guarded as ForeignArray
import qualified Foreign.Marshal.Utils as Marshal
import Foreign.C.Types (CInt, CChar)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Ptr (Ptr, nullPtr)
import Foreign.Storable (Storable)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)

import Data.Complex (Complex)
import Data.Tuple.HT (mapSnd)
import Data.Bool.HT (if')


values ::
   (Shape.C height, Shape.C width, Class.Floating a) =>
   General height width a -> Vector ZeroInt (RealOf a)
values =
   valuesGen $ \extent ->
      zeroInt $
      min
         (Shape.size $ Extent.height extent)
         (Shape.size $ Extent.width extent)

valuesTall ::
   (Extent.C vert, Shape.C height, Shape.C width, Class.Floating a) =>
   Full vert Extent.Small height width a -> Vector width (RealOf a)
valuesTall = valuesGen Extent.width

valuesWide ::
   (Extent.C horiz, Shape.C height, Shape.C width, Class.Floating a) =>
   Full Extent.Small horiz height width a -> Vector height (RealOf a)
valuesWide = valuesTall . Matrix.transpose

valuesGen ::
   (Extent.C vert, Extent.C horiz, Shape.C width, Shape.C height,
    Shape.C shape, Class.Floating a) =>
   (Extent vert horiz height width -> shape) ->
   Full vert horiz height width a -> Vector shape (RealOf a)
valuesGen resultShape =
   runTakeDiagonal $
   Class.switchFloating
      (TakeDiagonal $ valuesAux resultShape)
      (TakeDiagonal $ valuesAux resultShape)
      (TakeDiagonal $ valuesAux resultShape)
      (TakeDiagonal $ valuesAux resultShape)

valuesAux ::
   (Extent.C vert, Extent.C horiz, Shape.C width, Shape.C height,
    Shape.C shape, Class.Floating a, RealOf a ~ ar, Storable ar) =>
   (Extent vert horiz height width -> shape) ->
   Full vert horiz height width a -> Vector shape ar
valuesAux resultShape (Array shape@(MatrixShape.Full _order extent) a) =
   Array.unsafeCreateWithSize (resultShape extent) $ \mn sPtr -> do
   let (m,n) = MatrixShape.dimensions shape
   let lda = m
   evalContT $ do
      jobuPtr <- Call.char 'N'
      jobvtPtr <- Call.char 'N'
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      aPtr <- copyToTemp (m*n) a
      ldaPtr <- Call.leadingDim lda
      let uPtr = nullPtr
      let vtPtr = nullPtr
      lduPtr <- Call.leadingDim m
      ldvtPtr <- Call.leadingDim n
      liftIO $
         withInfo "gesvd" $ \infoPtr ->
         gesvd jobuPtr jobvtPtr mPtr nPtr
            aPtr ldaPtr sPtr uPtr lduPtr vtPtr ldvtPtr mn infoPtr


determinantAbsolute ::
   (Shape.C height, Shape.C width, Class.Floating a) =>
   General height width a -> RealOf a
determinantAbsolute =
   getDeterminant $
   Class.switchFloating
      (Determinant determinantAbsoluteAux)
      (Determinant determinantAbsoluteAux)
      (Determinant determinantAbsoluteAux)
      (Determinant determinantAbsoluteAux)

determinantAbsoluteAux ::
   (Shape.C height, Shape.C width,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General height width a -> ar
determinantAbsoluteAux =
   either (Vector.product . valuesTall) (const zero)
   .
   Matrix.caseTallWide


decompose ::
   (Shape.C height, Shape.C width, Class.Floating a) =>
   General height width a ->
   (Square height a, Vector ZeroInt (RealOf a), Square width a)
decompose =
   getDecompose $
   Class.switchFloating
      (Decompose decomposeAux)
      (Decompose decomposeAux)
      (Decompose decomposeAux)
      (Decompose decomposeAux)

newtype Decompose m f v g a =
   Decompose {getDecompose :: m a -> (f a, v (RealOf a), g a)}

decomposeAux ::
   (Shape.C height, Shape.C width,
    Class.Floating a, RealOf a ~ ar, Storable ar) =>
   General height width a ->
   (Square height a, Vector ZeroInt ar, Square width a)
decomposeAux arr@(Array shape@(MatrixShape.Full order extent) a) =

   let (height,width) = Extent.dimensions extent
       (m,n) = MatrixShape.dimensions shape
       mn = min m n

   in (if' (mn==0)
         (Square.identityFromHeight arr,
          Vector.autoFromList [],
          Square.identityFromWidth arr)) $
      (\(u,(s,vt)) -> (u,s,vt)) $
      Array.unsafeCreateWithSizeAndResult (MatrixShape.square order height) $
         \ _ uPtr0 ->
      ArrayIO.unsafeCreateWithSizeAndResult (zeroInt mn) $ \ _ sPtr ->
      ArrayIO.unsafeCreate (MatrixShape.square order width) $ \vtPtr0 ->

   evalContT $ do
      let (uPtr,vtPtr) = swapOnRowMajor order (uPtr0,vtPtr0)
      let lda = m
      jobuPtr <- Call.char 'A'
      jobvtPtr <- Call.char 'A'
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      aPtr <- copyToTemp (m*n) a
      ldaPtr <- Call.leadingDim lda
      lduPtr <- Call.leadingDim m
      ldvtPtr <- Call.leadingDim n
      liftIO $
         withInfo "gesvd" $ \infoPtr ->
         gesvd jobuPtr jobvtPtr mPtr nPtr
            aPtr ldaPtr sPtr uPtr lduPtr vtPtr ldvtPtr mn infoPtr


decomposeWide ::
   (Extent.C vert, Shape.C height, Shape.C width, Class.Floating a) =>
   Full Extent.Small vert height width a ->
   (Square height a, Vector height (RealOf a),
      Full Extent.Small vert height width a)
decomposeWide a =
   let (u,s,vt) = decomposeTall $ Matrix.transpose a
   in  (Square.transpose vt, s, Matrix.transpose u)

decomposeTall ::
   (Extent.C horiz, Shape.C height, Shape.C width, Class.Floating a) =>
   Full horiz Extent.Small height width a ->
   (Full horiz Extent.Small height width a,
      Vector width (RealOf a), Square width a)
decomposeTall =
   getDecompose $
   Class.switchFloating
      (Decompose decomposeThin)
      (Decompose decomposeThin)
      (Decompose decomposeThin)
      (Decompose decomposeThin)

decomposeThin ::
   (Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a, RealOf a ~ ar, Storable ar) =>
   Full horiz Extent.Small height width a ->
   (Full horiz Extent.Small height width a, Vector width ar, Square width a)
decomposeThin (Array (MatrixShape.Full order extent) a) =
   let (height,width) = Extent.dimensions extent
   in (\(u,(s,vt)) -> (u,s,vt)) $
      Array.unsafeCreateWithSizeAndResult (MatrixShape.Full order extent) $
         \ _ uPtr0 ->
      ArrayIO.unsafeCreateWithSizeAndResult width $ \ _ sPtr ->
      ArrayIO.unsafeCreate (MatrixShape.square order width) $ \vtPtr0 ->

   evalContT $ do
      let ((m,uPtr),(n,vtPtr)) =
            swapOnRowMajor order
               ((Shape.size height, uPtr0), (Shape.size width, vtPtr0))
      let mn = min m n
      let lda = m
      jobuPtr <- Call.char 'S'
      jobvtPtr <- Call.char 'S'
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      aPtr <- copyToTemp (m*n) a
      ldaPtr <- Call.leadingDim lda
      lduPtr <- Call.leadingDim m
      ldvtPtr <- Call.leadingDim mn
      liftIO $
         withInfo "gesvd" $ \infoPtr ->
         gesvd jobuPtr jobvtPtr mPtr nPtr
            aPtr ldaPtr sPtr uPtr lduPtr vtPtr ldvtPtr mn infoPtr


type GESVD_ ar a =
   Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt ->
   Ptr a -> Ptr CInt -> Ptr ar ->
   Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Int -> Ptr CInt -> IO ()

newtype GESVD a = GESVD {getGESVD :: GESVD_ (RealOf a) a}

gesvd :: Class.Floating a => GESVD_ (RealOf a) a
gesvd =
   getGESVD $
   Class.switchFloating
      (GESVD gesvdReal)
      (GESVD gesvdReal)
      (GESVD gesvdComplex)
      (GESVD gesvdComplex)

gesvdReal :: (Class.Real a) => GESVD_ a a
gesvdReal jobuPtr jobvtPtr mPtr nPtr
      aPtr ldaPtr sPtr uPtr lduPtr vtPtr ldvtPtr _mn infoPtr =
   withAutoWorkspace $ \workPtr lworkPtr ->
   LapackReal.gesvd jobuPtr jobvtPtr
      mPtr nPtr aPtr ldaPtr sPtr uPtr lduPtr vtPtr ldvtPtr
      workPtr lworkPtr infoPtr

gesvdComplex :: (Class.Real a) => GESVD_ a (Complex a)
gesvdComplex jobuPtr jobvtPtr
      mPtr nPtr aPtr ldaPtr sPtr uPtr lduPtr vtPtr ldvtPtr mn infoPtr =
   ForeignArray.alloca (5*mn) $ \rworkPtr ->
   withAutoWorkspace $ \workPtr lworkPtr ->
   LapackComplex.gesvd jobuPtr jobvtPtr
      mPtr nPtr aPtr ldaPtr sPtr uPtr lduPtr vtPtr ldvtPtr
      workPtr lworkPtr rworkPtr infoPtr


leastSquaresMinimumNormRCond ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Shape.C nrhs, Class.Floating a) =>
   RealOf a ->
   Full horiz vert height width a ->
   Full vert horiz height nrhs a ->
   (Int, Full vert horiz width nrhs a)
leastSquaresMinimumNormRCond rcond
      (Array (MatrixShape.Full orderA extentA) a)
      (Array (MatrixShape.Full orderB extentB) b) =
   case Extent.fuse (Extent.transpose extentA) extentB of
      Nothing -> error "leastSquaresMinimumNorm: height shapes mismatch"
      Just extent ->
         let widthA = Extent.width extentA
             (height,widthB) = Extent.dimensions extentB
             shapeX = MatrixShape.Full ColumnMajor extent
             m = Shape.size height
             n = Shape.size widthA
             nrhs = Shape.size widthB
         in if m == 0
               then (0, Vector.constant shapeX zero)
               else
                  if nrhs == 0
                     then
                        (fst $ unsafePerformIO $
                         case Vector.constant height zero of
                           Array _ b1 ->
                              leastSquaresMinimumNormIO rcond
                                 (MatrixShape.general ColumnMajor widthA ())
                                 orderA a orderB b1 m n 1,
                         Vector.constant shapeX zero)
                     else
                        unsafePerformIO $
                        leastSquaresMinimumNormIO rcond shapeX
                           orderA a orderB b m n nrhs

leastSquaresMinimumNormIO ::
   (Shape.C sh, Class.Floating a) =>
   RealOf a -> sh ->
   Order -> ForeignPtr a ->
   Order -> ForeignPtr a ->
   Int -> Int -> Int -> IO (Int, Array sh a)
leastSquaresMinimumNormIO rcond shapeX orderA a orderB b m n nrhs =
   createHigherArray shapeX m n nrhs $ \(tmpPtr,ldtmp) -> do

   let mn = min m n
   let aSize = m*n
   let lda = m
   evalContT $ do
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      nrhsPtr <- Call.cint nrhs
      aPtr <- Call.allocaArray aSize
      liftIO $ withForeignPtr a $ \asrcPtr ->
         copyToColumnMajor orderA m n asrcPtr aPtr
      ldaPtr <- Call.leadingDim lda
      ldtmpPtr <- Call.leadingDim ldtmp
      liftIO $ withForeignPtr b $ \bPtr ->
         copyToSubColumnMajor orderB m nrhs bPtr ldtmp tmpPtr

      rankPtr <- Call.alloca
      liftIO $
         withInfo "gelss" $ \infoPtr ->
         gelss mPtr nPtr nrhsPtr aPtr ldaPtr tmpPtr ldtmpPtr rcond
            rankPtr mn infoPtr

      liftIO $ peekCInt rankPtr


type GELSS_ ar a =
   Ptr CInt -> Ptr CInt -> Ptr CInt ->
   Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt ->
   ar -> Ptr CInt -> Int -> Ptr CInt -> IO ()

newtype GELSS a = GELSS {getGELSS :: GELSS_ (RealOf a) a}

gelss :: Class.Floating a => GELSS_ (RealOf a) a
gelss =
   getGELSS $
   Class.switchFloating
      (GELSS gelssReal)
      (GELSS gelssReal)
      (GELSS gelssComplex)
      (GELSS gelssComplex)

gelssReal :: (Class.Real a) => GELSS_ a a
gelssReal mPtr nPtr nrhsPtr aPtr ldaPtr bPtr ldbPtr rcond
      rankPtr mn infoPtr =
   Marshal.with rcond $ \rcondPtr ->
   ForeignArray.alloca mn $ \sPtr ->
   withAutoWorkspace $ \workPtr lworkPtr ->
   LapackReal.gelss
      mPtr nPtr nrhsPtr aPtr ldaPtr bPtr ldbPtr sPtr rcondPtr
      rankPtr workPtr lworkPtr infoPtr

gelssComplex :: (Class.Real a) => GELSS_ a (Complex a)
gelssComplex mPtr nPtr nrhsPtr aPtr ldaPtr bPtr ldbPtr rcond
      rankPtr mn infoPtr =
   Marshal.with rcond $ \rcondPtr ->
   ForeignArray.alloca mn $ \sPtr ->
   ForeignArray.alloca (5*mn) $ \rworkPtr ->
   withAutoWorkspace $ \workPtr lworkPtr ->
   LapackComplex.gelss
      mPtr nPtr nrhsPtr aPtr ldaPtr bPtr ldbPtr sPtr rcondPtr
      rankPtr workPtr lworkPtr rworkPtr infoPtr


pseudoInverseRCond ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Eq width, Class.Floating a) =>
   RealOf a ->
   Full vert horiz height width a ->
   (Int, Full horiz vert width height a)
pseudoInverseRCond =
   getPseudoInverseRCond $
   Class.switchFloating
      (PseudoInverseRCond pseudoInverseRCondAux)
      (PseudoInverseRCond pseudoInverseRCondAux)
      (PseudoInverseRCond pseudoInverseRCondAux)
      (PseudoInverseRCond pseudoInverseRCondAux)

newtype PseudoInverseRCond f g a =
   PseudoInverseRCond {
      getPseudoInverseRCond :: RealOf a -> f a -> (Int, g a)
   }

pseudoInverseRCondAux ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Eq width,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   ar ->
   Full vert horiz height width a ->
   (Int, Full horiz vert width height a)
pseudoInverseRCondAux rcond =
   getPseudoInverseExtent $
   Extent.switchTagPair
      (PseudoInverseExtent $ pseudoInverseRCondWide rcond)
      (PseudoInverseExtent $ pseudoInverseRCondWide rcond)
      (PseudoInverseExtent $ pseudoInverseRCondTall rcond)
      (PseudoInverseExtent $
         either
            (mapSnd Matrix.fromFull . pseudoInverseRCondTall rcond)
            (mapSnd Matrix.fromFull . pseudoInverseRCondWide rcond)
         .
         Matrix.caseTallWide)

newtype PseudoInverseExtent height width a vert horiz =
   PseudoInverseExtent {
      getPseudoInverseExtent ::
         Full vert horiz height width a ->
         (Int, Full horiz vert width height a)
   }

pseudoInverseRCondWide ::
   (Extent.C horiz, Shape.C height, Eq height, Shape.C width, Eq width,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   RealOf a ->
   Full Extent.Small horiz height width a ->
   (Int, Full horiz Extent.Small width height a)
pseudoInverseRCondWide rcond a =
   let (u,s,vt) = decomposeWide a
       (rank,recipS) = recipSigma rcond s
   in  (rank,
        Matrix.multiply (Matrix.adjoint vt) $
        scaleRowsReal recipS $ Square.toFull $ Square.adjoint u)

pseudoInverseRCondTall ::
   (Extent.C vert, Shape.C height, Eq height, Shape.C width, Eq width,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   RealOf a ->
   Full vert Extent.Small height width a ->
   (Int, Full Extent.Small vert width height a)
pseudoInverseRCondTall rcond a =
   let (u,s,vt) = decomposeTall a
       (rank,recipS) = recipSigma rcond s
   in  (rank,
        Matrix.multiply (Square.toFull $ Square.adjoint vt) $
        scaleRowsReal recipS $ Matrix.adjoint u)


recipSigma ::
   (Shape.C sh, Class.Real a) => a -> Array sh a -> (Int, Array sh a)
recipSigma rcond sigmas =
   case Array.toList sigmas of
      [] -> (0, sigmas)
      0:_ -> (0, sigmas)
      xs@(x:_) ->
         let smin = x * rcond
         in (length (takeWhile (>=smin) xs),
             Array.map (\s -> if s>=smin then recip s else 0) sigmas)


withInfo :: String -> (Ptr CInt -> IO ()) -> IO ()
withInfo = Private.withInfo "%d superdiagonals did not converge"
