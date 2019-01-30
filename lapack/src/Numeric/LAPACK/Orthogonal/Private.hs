{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Orthogonal.Private where

import qualified Numeric.LAPACK.Matrix.Private as Matrix
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as ExtentPriv
import qualified Numeric.LAPACK.Matrix.Extent as Extent
import qualified Numeric.LAPACK.Split as Split
import Numeric.LAPACK.Matrix.Triangular.Basic (Upper)
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor, ColumnMajor), sideSwapFromOrder)
import Numeric.LAPACK.Matrix.Extent.Private (Extent)
import Numeric.LAPACK.Matrix.Private
         (Full, ZeroInt, zeroInt,
          Transposition(NonTransposed, Transposed),
          Conjugation(NonConjugated, Conjugated),
          Inversion(NonInverted, Inverted), flipInversion)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Format (Format(format))
import Numeric.LAPACK.Scalar (RealOf, zero, isZero, absolute, conjugate)
import Numeric.LAPACK.Private
         (fill, copySubMatrix, copyBlock, conjugateToTemp,
          withAutoWorkspaceInfo, errorCodeMsg)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked.Monadic as ArrayIO
import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.Marshal.Array (advancePtr)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Monad (when)
import Control.Applicative (Const(Const,getConst), liftA3)

import qualified Data.List as List


data Householder vert horiz height width a =
   Householder {
      tau_ :: Vector ZeroInt a,
      split_ ::
         Array
            (MatrixShape.Split MatrixShape.Reflector vert horiz height width) a
   } deriving (Show)

type General = Householder Extent.Big Extent.Big
type Tall = Householder Extent.Big Extent.Small
type Wide = Householder Extent.Small Extent.Big
type Square sh = Householder Extent.Small Extent.Small sh sh


extent_ ::
   Householder vert horiz height width a ->
   Extent vert horiz height width
extent_ = MatrixShape.splitExtent . Array.shape . split_

mapExtent ::
   (Extent.C vertA, Extent.C horizA) =>
   (Extent.C vertB, Extent.C horizB) =>
   Extent.Map vertA horizA vertB horizB height width ->
   Householder vertA horizA height width a ->
   Householder vertB horizB height width a
mapExtent f (Householder tau split) =
   Householder tau $ Array.mapShape (MatrixShape.splitMapExtent f) split

caseTallWide ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
   Householder vert horiz height width a ->
   Either (Tall height width a) (Wide height width a)
caseTallWide (Householder tau (Array shape a)) =
   either
      (Left . Householder tau . flip Array a)
      (Right . Householder tau . flip Array a) $
   MatrixShape.caseTallWideSplit shape


instance
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
      Format (Householder vert horiz height width a) where
   format fmt (Householder tau m) = format fmt (tau, m)

fromMatrix ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Full vert horiz height width a ->
   Householder vert horiz height width a
fromMatrix (Array shape@(MatrixShape.Full order extent) a) =
   let (m,n) = MatrixShape.dimensions shape
   in uncurry Householder $
      Array.unsafeCreateWithSizeAndResult (zeroInt $ min m n) $ \_ tauPtr ->
      ArrayIO.unsafeCreate
         (MatrixShape.Split MatrixShape.Reflector order extent) $ \qrPtr ->

   evalContT $ do
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim m
      liftIO $ do
         copyBlock (m*n) aPtr qrPtr
         case order of
            RowMajor ->
               withAutoWorkspaceInfo errorCodeMsg "gelqf" $
                  LapackGen.gelqf mPtr nPtr qrPtr ldaPtr tauPtr
            ColumnMajor ->
               withAutoWorkspaceInfo errorCodeMsg "geqrf" $
                  LapackGen.geqrf mPtr nPtr qrPtr ldaPtr tauPtr

determinantR ::
   (Extent.C vert, Shape.C height, Shape.C width, Class.Floating a) =>
   Householder vert Extent.Small height width a -> a
determinantR = Split.determinantR . split_

{-
For complex numbers LAPACK uses not exactly reflections,
i.e. the determinants of the primitive transformations are not necessarily -1.

It holds: det(I-tau*v*v^H) = 1-tau*v^H*v
   because of https://en.wikipedia.org/wiki/Sylvester's_determinant_theorem
   simple proof from: https://en.wikipedia.org/wiki/Matrix_determinant_lemma
   I  0 . I+u*vt u .  I  0  =  I+u*vt     u      .  I  0 = I u
   vt 1     0    1   -vt 1     vt+vt*u*vt vt*u+1   -vt 1   0 vt*u+1

We already know:
   v^H*v is real and greater or equal to 1, because v[i] = 1,
   and determinant has absolute value 1.

Let k = v^H*v.
For which real k lies 1-tau*k on the unit circle?

   (1-taur*k)^2 + (taui*k)^2 = 1
   1-2*taur*k+(taur^2+taui^2)*k^2 = 1
   (taur^2 + taui^2)*k^2 - 2*taur*k = 0   (k/=0)
   (taur^2 + taui^2)*k - 2*taur = 0
   k = 2*taur / (taur^2 + taui^2)

   1-tau*k
      = (taur^2 + taui^2 - tau*2*taur) / (taur^2 + taui^2)
      = (taur^2 + taui^2 - 2*(taur+i*taui)*taur) / (taur^2 + taui^2)
      = (-taur^2 + taui^2 - 2*(i*taui)*taur) / (taur^2 + taui^2)
      = -(taur + i*taui)^2 / (taur^2 + taui^2)
-}
determinant ::
   (Shape.C sh, Class.Floating a) =>
   Square sh a -> a
determinant (Householder tau split) =
   List.foldl' (*) (Split.determinantR split) $
   (case MatrixShape.splitOrder $ Array.shape split of
      RowMajor -> map conjugate
      ColumnMajor -> id) $
   map (negate.(^(2::Int)).signum) $
   filter (not . isZero) $ Array.toList tau

determinantAbsolute ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Householder vert horiz height width a -> RealOf a
determinantAbsolute =
   absolute . either determinantR (const zero) . caseTallWide


leastSquares ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Eq width, Shape.C nrhs,
    Class.Floating a) =>
   Householder horiz Extent.Small height width a ->
   Full vert horiz height nrhs a ->
   Full vert horiz width nrhs a
leastSquares qr =
   tallSolveR NonTransposed NonConjugated qr . tallMultiplyQAdjoint qr

{- |
@
HH.minimumNorm (HH.fromMatrix a) b
==
Ortho.minimumNorm (adjoint a) b
@
-}
minimumNorm ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Eq width, Shape.C nrhs,
    Class.Floating a) =>
   Householder vert Extent.Small width height a ->
   Full vert horiz height nrhs a ->
   Full vert horiz width nrhs a
minimumNorm qr = tallMultiplyQ qr . tallSolveR Transposed Conjugated qr

-- cf. Matrix.takeRows
takeRows ::
   (Extent.C vert, Extent.C horiz,
    Eq fuse, Shape.C fuse, Shape.C height, Shape.C width, Class.Floating a) =>
   Extent Extent.Small horiz height fuse ->
   Full vert horiz fuse width a ->
   Full vert horiz height width a
takeRows extentA (Array (MatrixShape.Full order extentB) b) =
   case Extent.fuse (ExtentPriv.generalizeWide extentA) extentB of
      Nothing -> error "Householder.takeRows: heights mismatch"
      Just extentC ->
         Array.unsafeCreateWithSize (MatrixShape.Full order extentC) $
            \blockSize cPtr ->
         withForeignPtr b $ \bPtr ->
         case order of
            RowMajor -> copyBlock blockSize bPtr cPtr
            ColumnMajor ->
               let n  = Shape.size $ Extent.width  extentB
                   mb = Shape.size $ Extent.height extentB
                   mc = Shape.size $ Extent.height extentC
               in  copySubMatrix mc n mb bPtr mc cPtr

addRows ::
   (Extent.C vert, Extent.C horiz,
    Eq fuse, Shape.C fuse, Shape.C height, Shape.C width, Class.Floating a) =>
   Extent vert Extent.Small height fuse ->
   Full vert horiz fuse width a ->
   Full vert horiz height width a
addRows extentA (Array shapeB@(MatrixShape.Full order extentB) b) =
   case Extent.fuse (ExtentPriv.generalizeTall extentA) extentB of
      Nothing -> error "Householder.addRows: heights mismatch"
      Just extentC ->
         Array.unsafeCreateWithSize (MatrixShape.Full order extentC) $
            \cSize cPtr ->
         withForeignPtr b $ \bPtr ->
         case order of
            RowMajor -> do
               let bSize = Shape.size shapeB
               copyBlock bSize bPtr cPtr
               fill zero (cSize - bSize) (advancePtr cPtr bSize)
            ColumnMajor -> do
               let n  = Shape.size $ Extent.width  extentB
                   mb = Shape.size $ Extent.height extentB
                   mc = Shape.size $ Extent.height extentC
               copySubMatrix mb n mb bPtr mc cPtr
               evalContT $ do
                  uploPtr <- Call.char 'A'
                  mPtr <- Call.cint (mc-mb)
                  nPtr <- Call.cint n
                  ldcPtr <- Call.leadingDim mc
                  zPtr <- Call.number zero
                  liftIO $
                     LapackGen.laset uploPtr mPtr nPtr zPtr zPtr
                        (advancePtr cPtr mb) ldcPtr


extractQ ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Householder vert horiz height width a -> Matrix.Square height a
extractQ
   (Householder tau (Array (MatrixShape.Split _ order extent) qr)) =
      extractQAux tau (Extent.width extent) order
         (Extent.square $ Extent.height extent) qr

tallExtractQ ::
   (Extent.C vert, Shape.C height, Shape.C width, Class.Floating a) =>
   Householder vert Extent.Small height width a ->
   Full vert Extent.Small height width a
tallExtractQ
   (Householder tau (Array (MatrixShape.Split _ order extent) qr)) =
      extractQAux tau (Extent.width extent) order extent qr


extractQAux ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Shape.C widthQR,
    Class.Floating a) =>
   Vector ZeroInt a -> widthQR ->
   Order -> Extent vert horiz height width -> ForeignPtr a ->
   Full vert horiz height width a
extractQAux (Array widthTau tau) widthQR order extent qr =
   Array.unsafeCreate (MatrixShape.Full order extent) $ \qPtr -> do

   let (height,width) = Extent.dimensions extent
   let m = Shape.size height
   let n = Shape.size width
   let k = Shape.size widthTau
   evalContT $ do
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      qrPtr <- ContT $ withForeignPtr qr
      tauPtr <- ContT $ withForeignPtr tau
      case order of
         RowMajor -> do
            ldaPtr <- Call.leadingDim n
            liftIO $ do
               copySubMatrix k m (Shape.size widthQR) qrPtr n qPtr
               withAutoWorkspaceInfo errorCodeMsg "unglq" $
                  LapackGen.unglq nPtr mPtr kPtr qPtr ldaPtr tauPtr
         ColumnMajor -> do
            ldaPtr <- Call.leadingDim m
            liftIO $ do
               copyBlock (m*k) qrPtr qPtr
               withAutoWorkspaceInfo errorCodeMsg "ungqr" $
                  LapackGen.ungqr mPtr nPtr kPtr qPtr ldaPtr tauPtr


tallMultiplyQ ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Shape.C fuse, Eq fuse,
    Class.Floating a) =>
   Householder vert Extent.Small height fuse a ->
   Full vert horiz fuse width a ->
   Full vert horiz height width a
tallMultiplyQ qr = multiplyQ NonInverted qr . addRows (extent_ qr)

tallMultiplyQAdjoint ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Shape.C fuse, Eq fuse, Class.Floating a) =>
   Householder horiz Extent.Small fuse height a ->
   Full vert horiz fuse width a ->
   Full vert horiz height width a
tallMultiplyQAdjoint qr =
   takeRows (Extent.transpose $ extent_ qr) . multiplyQ Inverted qr


multiplyQ ::
   (Extent.C vertA, Extent.C horizA, Shape.C widthA,
    Extent.C vertB, Extent.C horizB, Shape.C widthB,
    Shape.C height, Eq height, Class.Floating a) =>
   Inversion ->
   Householder vertA horizA height widthA a ->
   Full vertB horizB height widthB a ->
   Full vertB horizB height widthB a
multiplyQ inverted
   (Householder
      (Array widthTau tau)
      (Array shapeA@(MatrixShape.Split _ orderA extentA) qr))
   (Array shapeB@(MatrixShape.Full orderB extentB) b) =

   Array.unsafeCreateWithSize shapeB $ \cSize cPtr -> do

   let (heightA,widthA) = Extent.dimensions extentA
   let (height,width) = Extent.dimensions extentB
   Call.assert "Householder.multiplyQ: height shapes mismatch"
      (heightA == height)

   let (side,(m,n)) =
         sideSwapFromOrder orderB (Shape.size height, Shape.size width)

   evalContT $ do
      sidePtr <- Call.char side
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      let k = Shape.size widthTau
      kPtr <- Call.cint k
      (transPtr,qrPtr,tauPtr) <-
         if orderA==orderB
           then
               liftA3 (,,)
                  (Call.char $ transposeFromInversion qr inverted)
                  (ContT $ withForeignPtr qr)
                  (ContT $ withForeignPtr tau)
           else
               liftA3 (,,)
                  (Call.char $
                   transposeFromInversion qr $ flipInversion inverted)
                  (conjugateToTemp (Shape.size shapeA) qr)
                  (conjugateToTemp k tau)
      bPtr <- ContT $ withForeignPtr b
      ldcPtr <- Call.leadingDim m
      liftIO $ copyBlock cSize bPtr cPtr
      case orderA of
         ColumnMajor -> do
            ldaPtr <- Call.leadingDim $ Shape.size heightA
            liftIO $ withAutoWorkspaceInfo errorCodeMsg "unmqr" $
               LapackGen.unmqr sidePtr transPtr
                  mPtr nPtr kPtr qrPtr ldaPtr tauPtr cPtr ldcPtr
         RowMajor -> do
            ldaPtr <- Call.leadingDim $ Shape.size widthA
            -- work-around for https://github.com/Reference-LAPACK/lapack/issues/260
            liftIO $ when (k>0) $
               withAutoWorkspaceInfo errorCodeMsg "unmlq" $
               LapackGen.unmlq sidePtr transPtr
                  mPtr nPtr kPtr qrPtr ldaPtr tauPtr cPtr ldcPtr

transposeFromInversion :: (Class.Floating a) => f a -> Inversion -> Char
transposeFromInversion qr Inverted = invChar qr
transposeFromInversion _ NonInverted = 'N'

invChar :: (Class.Floating a) => f a -> Char
invChar f = getConst $ asFuncTypeOf f inverseChar

asFuncTypeOf :: f a -> g a -> g a
asFuncTypeOf = const id

inverseChar :: (Class.Floating a) => Const Char a
inverseChar =
   Class.switchFloating (Const 'T') (Const 'T') (Const 'C') (Const 'C')


extractR ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Householder vert horiz height width a ->
   Full vert horiz height width a
extractR = Split.extractTriangle (Right MatrixShape.Triangle) . split_

tallExtractR ::
   (Extent.C vert, Shape.C height, Shape.C width, Class.Floating a) =>
   Householder vert Extent.Small height width a -> Upper width a
tallExtractR = Split.tallExtractR . split_

tallMultiplyR ::
   (Extent.C vertA, Extent.C vert, Extent.C horiz, Shape.C height, Eq height,
    Shape.C heightA, Shape.C widthB, Class.Floating a) =>
   Transposition ->
   Householder vertA Extent.Small heightA height a ->
   Full vert horiz height widthB a ->
   Full vert horiz height widthB a
tallMultiplyR transposed = Split.tallMultiplyR transposed . split_

tallSolveR ::
   (Extent.C vertA, Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Eq width, Shape.C nrhs, Class.Floating a) =>
   Transposition -> Conjugation ->
   Householder vertA Extent.Small height width a ->
   Full vert horiz width nrhs a -> Full vert horiz width nrhs a
tallSolveR transposed conjugated =
   Split.tallSolveR transposed conjugated . split_
