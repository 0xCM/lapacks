{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.BandedHermitianPositiveDefinite.Linear (
   solve,
   solveDecomposed,
   decompose,
   determinant,
   ) where

import qualified Numeric.LAPACK.Matrix.Banded.Basic as Banded
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Linear.Private (solver)
import Numeric.LAPACK.Matrix.BandedHermitian.Basic (BandedHermitian)
import Numeric.LAPACK.Matrix.Hermitian.Private (Determinant(..))
import Numeric.LAPACK.Matrix.Triangular.Private (copyTriangleToTemp)
import Numeric.LAPACK.Matrix.Shape.Private (uploFromOrder)
import Numeric.LAPACK.Matrix.Private (Full, Conjugation(Conjugated))
import Numeric.LAPACK.Scalar (RealOf, realPart)
import Numeric.LAPACK.Private (copyBlock, withInfo, rankMsg, definiteMsg)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num (integralFromProxy)
import Type.Base.Proxy (Proxy(Proxy))

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.ForeignPtr (withForeignPtr)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)


solve ::
   (Unary.Natural offDiag, Shape.C size, Eq size,
    Extent.C vert, Extent.C horiz, Shape.C nrhs, Class.Floating a) =>
   BandedHermitian offDiag size a ->
   Full vert horiz size nrhs a -> Full vert horiz size nrhs a
solve (Array (MatrixShape.BandedHermitian numOff orderA shA) a) =
   solver "BandedHermitian.solve" shA $ \n nPtr nrhsPtr xPtr ldxPtr -> do
      uploPtr <- Call.char $ uploFromOrder orderA
      let k = integralFromProxy numOff
      let lda = k+1
      kPtr <- Call.cint k
      aPtr <- copyTriangleToTemp Conjugated orderA (n*lda) a
      ldaPtr <- Call.leadingDim lda
      liftIO $
         withInfo definiteMsg "pbsv" $
            LapackGen.pbsv uploPtr nPtr kPtr nrhsPtr aPtr ldaPtr xPtr ldxPtr

{- |
> solve a b == solveDecomposed (decompose a) b
> solve (covariance u) b == solveDecomposed u b
-}
solveDecomposed ::
   (Unary.Natural offDiag, Shape.C size, Eq size,
    Extent.C vert, Extent.C horiz, Shape.C nrhs, Class.Floating a) =>
   Banded.Upper offDiag size a ->
   Full vert horiz size nrhs a -> Full vert horiz size nrhs a
solveDecomposed (Array (MatrixShape.Banded (_zero,numOff) orderA shA) a) =
   solver "BandedHermitian.solveDecomposed" (Extent.squareSize shA) $
         \n nPtr nrhsPtr xPtr ldxPtr -> do
      uploPtr <- Call.char $ uploFromOrder orderA
      let k = integralFromProxy numOff
      let lda = k+1
      kPtr <- Call.cint k
      aPtr <- copyTriangleToTemp Conjugated orderA (n*lda) a
      ldaPtr <- Call.leadingDim lda
      liftIO $
         withInfo rankMsg "pbtrs" $
            LapackGen.pbtrs uploPtr nPtr kPtr nrhsPtr aPtr ldaPtr xPtr ldxPtr


{- |
Cholesky decomposition
-}
decompose ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   BandedHermitian offDiag size a -> Banded.Upper offDiag size a
decompose (Array (MatrixShape.BandedHermitian numOff order sh) a) =
   Array.unsafeCreateWithSize
      (MatrixShape.bandedSquare (Proxy,numOff) order sh) $ \bSize bPtr -> do
   evalContT $ do
      let k = integralFromProxy numOff
      uploPtr <- Call.char $ uploFromOrder order
      nPtr <- Call.cint $ Shape.size sh
      kPtr <- Call.cint k
      aPtr <- ContT $ withForeignPtr a
      ldbPtr <- Call.leadingDim $ k+1
      liftIO $ do
         copyBlock bSize aPtr bPtr
         withInfo definiteMsg "pbtrf" $
            LapackGen.pbtrf uploPtr nPtr kPtr bPtr ldbPtr


determinant ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   BandedHermitian offDiag size a -> RealOf a
determinant =
   getDeterminant $
   Class.switchFloating
      (Determinant determinantAux) (Determinant determinantAux)
      (Determinant determinantAux) (Determinant determinantAux)

determinantAux ::
   (Unary.Natural offDiag, Shape.C size,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian offDiag size a -> ar
determinantAux =
   (^(2::Int)) . product . map realPart . Array.toList .
   Banded.takeDiagonal . decompose
