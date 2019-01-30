{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.HermitianPositiveDefinite.Linear (
   solve,
   solveDecomposed,
   inverse,
   decompose,
   determinant,
   ) where

import Numeric.LAPACK.Matrix.Hermitian.Basic (Hermitian)
import Numeric.LAPACK.Matrix.Triangular.Basic (Upper, takeDiagonal)
import Numeric.LAPACK.Matrix.Private (Full)

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Hermitian.Private (Determinant(..))
import Numeric.LAPACK.Matrix.Triangular.Private (copyTriangleToTemp)
import Numeric.LAPACK.Matrix.Shape.Private (NonUnit(NonUnit), uploFromOrder)
import Numeric.LAPACK.Matrix.Private (Conjugation(Conjugated))
import Numeric.LAPACK.Linear.Private (solver)
import Numeric.LAPACK.Scalar (RealOf, realPart)
import Numeric.LAPACK.Private (copyBlock, withInfo, rankMsg, definiteMsg)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))
import Data.Array.Comfort.Shape (triangleSize)

import Foreign.ForeignPtr (withForeignPtr)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)


solve ::
   (Extent.C vert, Extent.C horiz,
    Shape.C sh, Eq sh, Shape.C nrhs, Class.Floating a) =>
   Hermitian sh a ->
   Full vert horiz sh nrhs a -> Full vert horiz sh nrhs a
solve (Array (MatrixShape.Hermitian orderA shA) a) =
   solver "Hermitian.solve" shA $ \n nPtr nrhsPtr xPtr ldxPtr -> do
      uploPtr <- Call.char $ uploFromOrder orderA
      apPtr <- copyTriangleToTemp Conjugated orderA (triangleSize n) a
      liftIO $
         withInfo definiteMsg "ppsv" $
            LapackGen.ppsv uploPtr nPtr nrhsPtr apPtr xPtr ldxPtr

{- |
> solve a b == solveDecomposed (decompose a) b
> solve (covariance u) b == solveDecomposed u b
-}
solveDecomposed ::
   (Extent.C vert, Extent.C horiz,
    Shape.C sh, Eq sh, Shape.C nrhs, Class.Floating a) =>
   Upper sh a -> Full vert horiz sh nrhs a -> Full vert horiz sh nrhs a
solveDecomposed (Array (MatrixShape.Triangular NonUnit _uplo orderA shA) a) =
   solver "Hermitian.solveDecomposed" shA $ \n nPtr nrhsPtr xPtr ldxPtr -> do
      uploPtr <- Call.char $ uploFromOrder orderA
      apPtr <- copyTriangleToTemp Conjugated orderA (triangleSize n) a
      liftIO $
         withInfo rankMsg "pptrs" $
            LapackGen.pptrs uploPtr nPtr nrhsPtr apPtr xPtr ldxPtr


inverse ::
   (Shape.C sh, Class.Floating a) => Hermitian sh a -> Hermitian sh a
inverse
   (Array shape@(MatrixShape.Hermitian order sh) a) =
      Array.unsafeCreateWithSize shape $ \triSize bPtr -> do
   evalContT $ do
      uploPtr <- Call.char $ uploFromOrder order
      nPtr <- Call.cint $ Shape.size sh
      aPtr <- ContT $ withForeignPtr a
      liftIO $ do
         copyBlock triSize aPtr bPtr
         withInfo definiteMsg "pptrf" $ LapackGen.pptrf uploPtr nPtr bPtr
         withInfo rankMsg "pptri" $ LapackGen.pptri uploPtr nPtr bPtr

{- |
Cholesky decomposition
-}
decompose ::
   (Shape.C sh, Class.Floating a) => Hermitian sh a -> Upper sh a
decompose
   (Array (MatrixShape.Hermitian order sh) a) =
      Array.unsafeCreateWithSize
         (MatrixShape.Triangular NonUnit MatrixShape.upper order sh) $
            \triSize bPtr -> do
   evalContT $ do
      uploPtr <- Call.char $ uploFromOrder order
      nPtr <- Call.cint $ Shape.size sh
      aPtr <- ContT $ withForeignPtr a
      liftIO $ do
         copyBlock triSize aPtr bPtr
         withInfo definiteMsg "pptrf" $ LapackGen.pptrf uploPtr nPtr bPtr


determinant ::
   (Shape.C sh, Class.Floating a) => Hermitian sh a -> RealOf a
determinant =
   getDeterminant $
   Class.switchFloating
      (Determinant determinantAux) (Determinant determinantAux)
      (Determinant determinantAux) (Determinant determinantAux)

determinantAux ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian sh a -> ar
determinantAux =
   (^(2::Int)) . product . map realPart . Array.toList . takeDiagonal . decompose
