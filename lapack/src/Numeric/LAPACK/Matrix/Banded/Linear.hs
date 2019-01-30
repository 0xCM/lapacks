{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.Banded.Linear (
   solve,
   solveColumnMajor,
   determinant,
   ) where

import qualified Numeric.LAPACK.Matrix.Banded.Basic as Banded
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Split as Split
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Linear.Private (solver, withDeterminantInfo, withInfo)
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor,ColumnMajor), transposeFromOrder)
import Numeric.LAPACK.Matrix.Private (Full)
import Numeric.LAPACK.Private (copySubMatrix)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num (integralFromProxy)

import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import System.IO.Unsafe (unsafePerformIO)

import Foreign.Marshal.Array (peekArray, advancePtr)
import Foreign.ForeignPtr (withForeignPtr)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)


solve ::
   (Unary.Natural sub, Unary.Natural super, Extent.C vert, Extent.C horiz,
    Shape.C sh, Eq sh, Shape.C nrhs, Class.Floating a) =>
   Banded.Square sub super sh a ->
   Full vert horiz sh nrhs a -> Full vert horiz sh nrhs a
solve (Array (MatrixShape.Banded numOff order extent) a) =
   solver "Banded.solve" (Extent.squareSize extent) $
         \n nPtr nrhsPtr xPtr ldxPtr -> do
      let (kl,ku) = MatrixShape.numOffDiagonals order numOff
      let k = kl+1+ku
      let ldab = kl+k
      transPtr <- Call.char $ transposeFromOrder order
      klPtr <- Call.cint kl
      kuPtr <- Call.cint ku
      aPtr <- ContT $ withForeignPtr a
      abPtr <- Call.allocaArray (n*ldab)
      ldabPtr <- Call.leadingDim ldab
      ipivPtr <- Call.allocaArray n
      liftIO $ do
         copySubMatrix k n k aPtr ldab (advancePtr abPtr kl)
         withInfo "gbtrf" $
            LapackGen.gbtrf nPtr nPtr klPtr kuPtr abPtr ldabPtr ipivPtr
         withInfo "gbtrs" $
            LapackGen.gbtrs transPtr nPtr klPtr kuPtr nrhsPtr
               abPtr ldabPtr ipivPtr xPtr ldxPtr

solveColumnMajor ::
   (Unary.Natural sub, Unary.Natural super, Extent.C vert, Extent.C horiz,
    Shape.C sh, Eq sh, Shape.C nrhs, Class.Floating a) =>
   Banded.Square sub super sh a ->
   Full vert horiz sh nrhs a -> Full vert horiz sh nrhs a
solveColumnMajor
      (Array (MatrixShape.Banded (sub,super) ColumnMajor extent) a) =
   solver "Banded.solve" (Extent.squareSize extent) $
         \n nPtr nrhsPtr xPtr ldxPtr -> do
      let kl = integralFromProxy sub
      let ku = integralFromProxy super
      let k = kl+1+ku
      let ldab = kl+k
      klPtr <- Call.cint kl
      kuPtr <- Call.cint ku
      aPtr <- ContT $ withForeignPtr a
      abPtr <- Call.allocaArray (n*ldab)
      ldabPtr <- Call.leadingDim ldab
      ipivPtr <- Call.allocaArray n
      liftIO $ do
         copySubMatrix k n k aPtr ldab (advancePtr abPtr kl)
         withInfo "gbsv" $
            LapackGen.gbsv nPtr klPtr kuPtr nrhsPtr
               abPtr ldabPtr ipivPtr xPtr ldxPtr
solveColumnMajor (Array (MatrixShape.Banded _ RowMajor _) _) =
   error "Linear.Banded.solveColumnMajor: RowMajor intentionally unimplemented"

determinant ::
   (Unary.Natural sub, Unary.Natural super, Shape.C sh, Class.Floating a) =>
   Banded.Square sub super sh a -> a
determinant (Array (MatrixShape.Banded numOff order extent) a) =
      unsafePerformIO $ do
   let n = Shape.size $ Extent.squareSize extent
   evalContT $ do
      let (kl,ku) = MatrixShape.numOffDiagonals order numOff
      let k = kl+1+ku
      let ldab = kl+k
      nPtr <- Call.cint n
      klPtr <- Call.cint kl
      kuPtr <- Call.cint ku
      aPtr <- ContT $ withForeignPtr a
      abPtr <- Call.allocaArray (n*ldab)
      ldabPtr <- Call.leadingDim ldab
      ipivPtr <- Call.allocaArray n
      liftIO $ do
         copySubMatrix k n k aPtr ldab (advancePtr abPtr kl)
         withDeterminantInfo "gbtrf"
            (LapackGen.gbtrf nPtr nPtr klPtr kuPtr abPtr ldabPtr ipivPtr)
            (do
               det <- Private.product n (advancePtr abPtr (kl+ku)) ldab
               ipiv <- peekArray n ipivPtr
               return $ if Split.oddPermutation ipiv then -det else det)
