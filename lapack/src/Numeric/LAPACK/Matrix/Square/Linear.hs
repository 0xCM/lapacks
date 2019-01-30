{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.Square.Linear (
   solve,
   inverse,
   determinant,
   ) where

import Numeric.LAPACK.Matrix.Private (Full, Square, argSquare)

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Split as Split
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Linear.Private
         (solver, withDeterminantInfo, withInfo, diagonalMsg)
import Numeric.LAPACK.Matrix.Shape.Private (transposeFromOrder)
import Numeric.LAPACK.Private
         (withAutoWorkspaceInfo, copyBlock, copyToTemp, copyToColumnMajor)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import System.IO.Unsafe (unsafePerformIO)

import Foreign.Marshal.Array (peekArray)
import Foreign.ForeignPtr (withForeignPtr)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Monad (when)


solve, _solve ::
   (Extent.C vert, Extent.C horiz,
    Shape.C sh, Eq sh, Shape.C nrhs, Class.Floating a) =>
   Square sh a -> Full vert horiz sh nrhs a -> Full vert horiz sh nrhs a
solve =
   argSquare $ \orderA shA a ->
   solver "Square.solve" shA $ \n nPtr nrhsPtr xPtr ldxPtr -> do
      transPtr <- Call.char $ transposeFromOrder orderA
      aPtr <- copyToTemp (n*n) a
      ldaPtr <- Call.leadingDim n
      ipivPtr <- Call.allocaArray n
      liftIO $ do
         withInfo "getrf" $
            LapackGen.getrf nPtr nPtr aPtr ldaPtr ipivPtr
         withInfo "getrs" $
            LapackGen.getrs transPtr nPtr nrhsPtr
               aPtr ldaPtr ipivPtr xPtr ldxPtr

_solve =
   argSquare $ \orderA shA a ->
   solver "Square.solve" shA $ \n nPtr nrhsPtr xPtr ldxPtr -> do
      aPtr <- ContT $ withForeignPtr a
      atmpPtr <- Call.allocaArray (n*n)
      ldaPtr <- Call.leadingDim n
      ipivPtr <- Call.allocaArray n
      liftIO $ do
         copyToColumnMajor orderA n n aPtr atmpPtr
         withInfo "gesv" $
            LapackGen.gesv nPtr nrhsPtr atmpPtr ldaPtr ipivPtr xPtr ldxPtr


inverse :: (Shape.C sh, Class.Floating a) => Square sh a -> Square sh a
inverse (Array shape@(MatrixShape.Full _order extent) a) =
      Array.unsafeCreateWithSize shape $ \blockSize bPtr -> do
   let n = Shape.size $ Extent.squareSize extent
   evalContT $ do
      nPtr <- Call.cint n
      aPtr <- ContT $ withForeignPtr a
      ldbPtr <- Call.leadingDim n
      ipivPtr <- Call.allocaArray n
      liftIO $ when (n>0) $ do
         copyBlock blockSize aPtr bPtr
         withInfo "getrf" $ LapackGen.getrf nPtr nPtr bPtr ldbPtr ipivPtr
         withAutoWorkspaceInfo diagonalMsg "getri" $
            LapackGen.getri nPtr bPtr ldbPtr ipivPtr


determinant :: (Shape.C sh, Class.Floating a) => Square sh a -> a
determinant = argSquare $ \_order sh a -> unsafePerformIO $ do
   let n = Shape.size sh
   evalContT $ do
      nPtr <- Call.cint n
      aPtr <- copyToTemp (n*n) a
      ldaPtr <- Call.leadingDim n
      ipivPtr <- Call.allocaArray n
      liftIO $ withDeterminantInfo "getrf"
         (LapackGen.getrf nPtr nPtr aPtr ldaPtr ipivPtr)
         (do
            det <- Private.product n aPtr (n+1)
            ipiv <- peekArray n ipivPtr
            return $ if Split.oddPermutation ipiv then -det else det)
