{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
module Numeric.LAPACK.Matrix.Triangular.Linear (
   solve,
   inverse,
   inverseGeneric,
   determinant,
   ) where

import qualified Numeric.LAPACK.Matrix.Banded.Linear as BandedLin
import qualified Numeric.LAPACK.Matrix.Banded.Basic as Banded
import qualified Numeric.LAPACK.Matrix.Symmetric.Private as Symmetric
import qualified Numeric.LAPACK.Matrix.Triangular.Private as Tri
import Numeric.LAPACK.Linear.Private (solver, withInfo)
import Numeric.LAPACK.Matrix.Triangular.Basic
         (Triangular, Symmetric, PowerDiag, takeDiagonal, strictNonUnitDiagonal)
import Numeric.LAPACK.Matrix.Private (Full)

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Shape.Private
         (transposeFromOrder, uploFromOrder, uploOrder, charFromTriDiag)
import Numeric.LAPACK.Matrix.Private (Conjugation(NonConjugated))
import Numeric.LAPACK.Private (copyBlock, copyToTemp)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))
import Data.Array.Comfort.Shape (triangleSize)

import System.IO.Unsafe (unsafePerformIO)

import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Ptr (Ptr)
import Foreign.Storable (peek)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)


solve ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Extent.C vert, Extent.C horiz,
    Shape.C sh, Eq sh, Shape.C nrhs, Class.Floating a) =>
   Triangular lo diag up sh a ->
   Full vert horiz sh nrhs a -> Full vert horiz sh nrhs a
solve =
   Tri.getMultiplyRight $
   MatrixShape.switchDiagUpLoSym
      (Tri.MultiplyRight $
       Tri.multiplyDiagonal
         "solve.diagonal: sizes mismatch"
         (MatrixShape.fullHeight . Array.shape)
         (BandedLin.solve . Banded.diagonal . takeDiagonal))
      (Tri.MultiplyRight solveTriangular)
      (Tri.MultiplyRight solveTriangular)
      (Tri.MultiplyRight $ solveSymmetric . strictNonUnitDiagonal)

solveTriangular ::
   (MatrixShape.UpLo lo up, MatrixShape.TriDiag diag,
    Extent.C vert, Extent.C horiz,
    Shape.C sh, Eq sh, Shape.C nrhs, Class.Floating a) =>
   Triangular lo diag up sh a ->
   Full vert horiz sh nrhs a -> Full vert horiz sh nrhs a
solveTriangular (Array (MatrixShape.Triangular diag uplo orderA shA) a) =
   solver "Triangular.solve" shA $ \n nPtr nrhsPtr xPtr ldxPtr -> do
      uploPtr <- Call.char $ uploFromOrder $ uploOrder uplo orderA
      transPtr <- Call.char $ transposeFromOrder orderA
      diagPtr <- Call.char $ charFromTriDiag diag
      apPtr <- copyToTemp (triangleSize n) a
      liftIO $
         withInfo "tptrs" $
            LapackGen.tptrs uploPtr transPtr diagPtr
               nPtr nrhsPtr apPtr xPtr ldxPtr

solveSymmetric ::
   (Extent.C vert, Extent.C horiz,
    Shape.C sh, Eq sh, Shape.C nrhs, Class.Floating a) =>
   Symmetric sh a ->
   Full vert horiz sh nrhs a -> Full vert horiz sh nrhs a
solveSymmetric (Array (MatrixShape.Triangular _diag _uplo orderA shA) a) =
   Symmetric.solve "Symmetric.solve" NonConjugated orderA shA a


inverse ::
   (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
    Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Triangular lo diag up sh a
inverse =
   Tri.getMap $
   MatrixShape.switchDiagUpLo
      (Tri.Map inverseDiagonal)
      (Tri.Map inverseTriangular)
      (Tri.Map inverseTriangular)

inverseGeneric ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a ->
   Triangular lo (PowerDiag lo up diag) up sh a
inverseGeneric =
   Tri.getPower $
   MatrixShape.switchDiagUpLoSym
      (Tri.Power inverseDiagonal)
      (Tri.Power inverseTriangular)
      (Tri.Power inverseTriangular)
      (Tri.Power $ inverseSymmetric . strictNonUnitDiagonal)

inverseDiagonal ::
   (MatrixShape.TriDiag diag, Shape.C sh, Class.Floating a) =>
   Tri.FlexDiagonal diag sh a -> Tri.FlexDiagonal diag sh a
inverseDiagonal = Tri.caseTriDiagArray id (Array.map recip)

inverseTriangular ::
   (MatrixShape.UpLo lo up, MatrixShape.TriDiag diag,
    Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Triangular lo diag up sh a
inverseTriangular (Array shape@(MatrixShape.Triangular diag uplo order sh) a) =
      Array.unsafeCreateWithSize shape $ \triSize bPtr ->
   evalContT $ do
      uploPtr <- Call.char $ uploFromOrder $ uploOrder uplo order
      diagPtr <- Call.char $ charFromTriDiag diag
      nPtr <- Call.cint $ Shape.size sh
      aPtr <- ContT $ withForeignPtr a
      liftIO $ do
         copyBlock triSize aPtr bPtr
         withInfo "tptri" $ LapackGen.tptri uploPtr diagPtr nPtr bPtr

inverseSymmetric ::
   (Shape.C sh, Class.Floating a) => Symmetric sh a -> Symmetric sh a
inverseSymmetric (Array shape@(MatrixShape.Triangular _diag _uplo order sh) a) =
   Array.unsafeCreateWithSize shape $
      Symmetric.inverse NonConjugated order (Shape.size sh) a


determinant ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a -> a
determinant =
   Tri.getMultiplyRight $
   MatrixShape.switchDiagUpLoSym
      (Tri.MultiplyRight determinantTriangular)
      (Tri.MultiplyRight determinantTriangular)
      (Tri.MultiplyRight determinantTriangular)
      (Tri.MultiplyRight $ determinantSymmetric . strictNonUnitDiagonal)

determinantTriangular ::
   (MatrixShape.DiagUpLo lo up, Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a -> a
determinantTriangular = product . Array.toList . takeDiagonal

determinantSymmetric ::
   (Shape.C sh, Class.Floating a) => Symmetric sh a -> a
determinantSymmetric (Array (MatrixShape.Triangular _diag _uplo order sh) a) =
   unsafePerformIO $
      Symmetric.determinant NonConjugated
         peekBlockDeterminant order (Shape.size sh) a

peekBlockDeterminant ::
   (Class.Floating a) => (Ptr a, Maybe (Ptr a, Ptr a)) -> IO a
peekBlockDeterminant (a0Ptr,ext) = do
   a0 <- peek a0Ptr
   case ext of
      Nothing -> return a0
      Just (a1Ptr,bPtr) -> do
         a1 <- peek a1Ptr
         b <- peek bPtr
         return (a0*a1 - b*b)
