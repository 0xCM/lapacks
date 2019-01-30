{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Orthogonal (
   leastSquares,
   minimumNorm,
   leastSquaresMinimumNormRCond,
   pseudoInverseRCond,

   determinant,
   determinantAbsolute,
   complement,

   householder,
   ) where

import qualified Numeric.LAPACK.Orthogonal.Private as HH

import qualified Numeric.LAPACK.Matrix.Square.Basic as Square
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Matrix.Extent.Kind as EK
import Numeric.LAPACK.Matrix.Square.Basic (Square)
import Numeric.LAPACK.Matrix.Private (Full, Tall, ZeroInt, zeroInt)
import Numeric.LAPACK.Matrix (transpose, dropColumns)

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix.Shape.Private (Order(RowMajor,ColumnMajor))
import Numeric.LAPACK.Scalar (RealOf, zero, absolute)
import Numeric.LAPACK.Private
         (lacgv, peekCInt,
          copySubMatrix, copyToTemp, copyToColumnMajor, copyToSubColumnMajor,
          withAutoWorkspaceInfo, rankMsg, errorCodeMsg, createHigherArray)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.LAPACK.FFI.Complex as LapackComplex
import qualified Numeric.LAPACK.FFI.Real as LapackReal
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import System.IO.Unsafe (unsafePerformIO)

import Foreign.Marshal.Array (pokeArray)
import Foreign.C.Types (CInt, CChar)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Ptr (Ptr)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)

import Data.Complex (Complex)
import Data.Tuple.HT (mapSnd)


{- |
If @x = leastSquares a b@
then @x@ minimizes @Vector.norm2 (multiply a x `sub` b)@.

Precondition: @a@ must have full rank and @height a >= width a@.
-}
leastSquares ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Shape.C nrhs, Class.Floating a) =>
   Full horiz Extent.Small height width a ->
   Full vert horiz height nrhs a ->
   Full vert horiz width nrhs a
leastSquares
   (Array shapeA@(MatrixShape.Full orderA extentA) a)
   (Array shapeB@(MatrixShape.Full orderB extentB) b) =

 case Extent.fuse (Extent.generalizeWide $ Extent.transpose extentA) extentB of
  Nothing -> error "leastSquares: height shapes mismatch"
  Just extent ->
      Array.unsafeCreate (MatrixShape.Full ColumnMajor extent) $ \xPtr -> do

   let widthA = Extent.width extentA
   let (height,widthB) = Extent.dimensions extentB
   let (m,n) = MatrixShape.dimensions shapeA
   let lda = m
   let nrhs = Shape.size widthB
   let ldb = Shape.size height
   let ldx = Shape.size widthA
   evalContT $ do
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      nrhsPtr <- Call.cint nrhs
      (transPtr,aPtr) <- transposeA orderA (m*n) a
      ldaPtr <- Call.leadingDim lda
      bPtr <- ContT $ withForeignPtr b
      ldbPtr <- Call.leadingDim ldb
      let bSize = Shape.size shapeB
      btmpPtr <- Call.allocaArray bSize
      liftIO $ copyToColumnMajor orderB ldb nrhs bPtr btmpPtr
      liftIO $ withAutoWorkspaceInfo rankMsg "gels" $
         LapackGen.gels transPtr
            mPtr nPtr nrhsPtr aPtr ldaPtr btmpPtr ldbPtr
      liftIO $ copySubMatrix ldx nrhs ldb btmpPtr ldx xPtr

{- |
The vector @x@ with @x = minimumNorm a b@
is the vector with minimal @Vector.norm2 x@
that satisfies @multiply a x == b@.

Precondition: @a@ must have full rank and @height a <= width a@.
-}
minimumNorm ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Shape.C nrhs, Class.Floating a) =>
   Full Extent.Small vert height width a ->
   Full vert horiz height nrhs a ->
   Full vert horiz width nrhs a
minimumNorm
   (Array shapeA@(MatrixShape.Full orderA extentA) a)
   (Array        (MatrixShape.Full orderB extentB) b) =

 case Extent.fuse (Extent.generalizeTall $ Extent.transpose extentA) extentB of
  Nothing -> error "minimumNorm: height shapes mismatch"
  Just extent ->
      Array.unsafeCreate (MatrixShape.Full ColumnMajor extent) $ \xPtr -> do

   let widthA = Extent.width extentA
   let (height,widthB) = Extent.dimensions extentB
   let (m,n) = MatrixShape.dimensions shapeA
   let lda = m
   let nrhs = Shape.size widthB
   let ldb = Shape.size height
   let ldx = Shape.size widthA
   evalContT $ do
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      nrhsPtr <- Call.cint nrhs
      (transPtr,aPtr) <- transposeA orderA (m*n) a
      ldaPtr <- Call.leadingDim lda
      bPtr <- ContT $ withForeignPtr b
      ldxPtr <- Call.leadingDim ldx
      liftIO $ copyToSubColumnMajor orderB ldb nrhs bPtr ldx xPtr
      liftIO $ withAutoWorkspaceInfo rankMsg "gels" $
         LapackGen.gels transPtr
            mPtr nPtr nrhsPtr aPtr ldaPtr xPtr ldxPtr


transposeA ::
   Class.Floating a =>
   Order -> Int -> ForeignPtr a -> ContT r IO (Ptr CChar, Ptr a)
transposeA order size a = do
   aPtr <- copyToTemp size a
   trans <-
      case order of
         RowMajor -> do
            sizePtr <- Call.cint size
            incPtr <- Call.cint 1
            liftIO $ lacgv sizePtr aPtr incPtr
            return $ HH.invChar a
         ColumnMajor -> return 'N'
   transPtr <- Call.char trans
   return (transPtr, aPtr)


{- |
If @x = leastSquaresMinimumNormRCond rcond a b@
then @x@ is the vector with minimum @Vector.norm2 x@
that minimizes @Vector.norm2 (multiply a x `sub` b)@.

Matrix @a@ can have any rank
but you must specify the reciprocal condition of the rank-truncated matrix.
-}
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
      Nothing -> error "leastSquaresMinimumNormRCond: height shapes mismatch"
      Just extent ->
         let widthA = Extent.width extentA
             (height,widthB) = Extent.dimensions extentB
             shapeX = MatrixShape.Full ColumnMajor extent
             m = Shape.size height
             n = Shape.size widthA
             nrhs = Shape.size widthB
         in  if m == 0
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

   let aSize = m*n
   let lda = m
   evalContT $ do
      aPtr <- ContT $ withForeignPtr a
      atmpPtr <- Call.allocaArray aSize
      liftIO $ copyToColumnMajor orderA m n aPtr atmpPtr
      ldaPtr <- Call.leadingDim lda
      ldtmpPtr <- Call.leadingDim ldtmp
      bPtr <- ContT $ withForeignPtr b
      liftIO $ copyToSubColumnMajor orderB m nrhs bPtr ldtmp tmpPtr
      jpvtPtr <- Call.allocaArray n
      liftIO $ pokeArray jpvtPtr (replicate n 0)
      rankPtr <- Call.alloca
      gelsy m n nrhs atmpPtr ldaPtr tmpPtr ldtmpPtr jpvtPtr rcond rankPtr
      liftIO $ peekCInt rankPtr


type GELSY_ r ar a =
   Int -> Int -> Int -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt ->
   Ptr CInt -> ar -> Ptr CInt -> ContT r IO ()

newtype GELSY r a = GELSY {getGELSY :: GELSY_ r (RealOf a) a}

gelsy :: (Class.Floating a) => GELSY_ r (RealOf a) a
gelsy =
   getGELSY $
   Class.switchFloating
      (GELSY gelsyReal)
      (GELSY gelsyReal)
      (GELSY gelsyComplex)
      (GELSY gelsyComplex)

gelsyReal :: (Class.Real a) => GELSY_ r a a
gelsyReal m n nrhs aPtr ldaPtr bPtr ldbPtr jpvtPtr rcond rankPtr = do
   mPtr <- Call.cint m
   nPtr <- Call.cint n
   nrhsPtr <- Call.cint nrhs
   rcondPtr <- Call.real rcond
   liftIO $ withAutoWorkspaceInfo errorCodeMsg "gelsy" $
      LapackReal.gelsy mPtr nPtr nrhsPtr
         aPtr ldaPtr bPtr ldbPtr jpvtPtr rcondPtr rankPtr

gelsyComplex :: (Class.Real a) => GELSY_ r a (Complex a)
gelsyComplex m n nrhs aPtr ldaPtr bPtr ldbPtr jpvtPtr rcond rankPtr = do
   mPtr <- Call.cint m
   nPtr <- Call.cint n
   nrhsPtr <- Call.cint nrhs
   rcondPtr <- Call.real rcond
   rworkPtr <- Call.allocaArray (2*n)
   liftIO $
      withAutoWorkspaceInfo errorCodeMsg "gelsy" $ \workPtr lworkPtr infoPtr ->
      LapackComplex.gelsy mPtr nPtr nrhsPtr
         aPtr ldaPtr bPtr ldbPtr jpvtPtr rcondPtr rankPtr
         workPtr lworkPtr rworkPtr infoPtr


pseudoInverseRCond ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Eq width, Class.Floating a) =>
   RealOf a ->
   Full vert horiz height width a ->
   (Int, Full horiz vert width height a)
pseudoInverseRCond rcond a =
   case Matrix.caseTallWide a of
      Left _ ->
         mapSnd transpose $
         leastSquaresMinimumNormRCond rcond (transpose a) $
         Square.toFull $ Square.identity $
         MatrixShape.fullWidth $ Array.shape a
      Right _ ->
         leastSquaresMinimumNormRCond rcond a $
         Square.toFull $ Square.identity $
         MatrixShape.fullHeight $ Array.shape a


{-
@(q,r) = householder a@
means that @q@ is unitary and @r@ is upper triangular and @a = multiply q r@.
-}
householder ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Full vert horiz height width a ->
   (Square height a, Full vert horiz height width a)
householder a =
   let hh = HH.fromMatrix a
   in  (HH.extractQ hh, HH.extractR hh)


determinant :: (Shape.C sh, Class.Floating a) => Square sh a -> a
determinant = HH.determinant . HH.fromMatrix

{-|
Gramian determinant -
works also for non-square matrices, but is sensitive to transposition.

> determinantAbsolute a = sqrt (Herm.determinant (Herm.covariance a))
-}
determinantAbsolute ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Full vert horiz height width a -> RealOf a
determinantAbsolute =
   absolute .
   either (HH.determinantR . HH.fromMatrix) (const zero) .
   Matrix.caseTallWide


{- |
For an m-by-n-matrix @a@ with m>=n
this function computes an m-by-(m-n)-matrix @b@
such that @Matrix.multiply (adjoint b) a@ is a zero matrix.
The function does not try to compensate a rank deficiency of @a@.
That is, @a|||b@ has full rank if and only if @a@ has full rank.

For full-rank matrices you might also call this @kernel@ or @nullspace@.
-}
complement ::
   (Shape.C height, Shape.C width, Class.Floating a) =>
   Tall height width a -> Tall height ZeroInt a
complement a =
   dropColumns (Shape.size $ MatrixShape.fullWidth $ Array.shape a) $
   Array.mapShape zeroIntWidth $ Square.toFull $
   HH.extractQ $ HH.fromMatrix $ Matrix.fromFull a

zeroIntWidth ::
   (Shape.C width) =>
   MatrixShape.Tall height width -> MatrixShape.Tall height ZeroInt
zeroIntWidth (MatrixShape.Full order (Extent.Extent o (EK.Tall height width))) =
   MatrixShape.Full order
      (Extent.Extent o (EK.Tall height (zeroInt $ Shape.size width)))
