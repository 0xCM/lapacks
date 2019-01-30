module Numeric.LAPACK.Linear.Private where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Matrix.Shape.Private (Order(ColumnMajor))
import Numeric.LAPACK.Matrix.Private (Full)
import Numeric.LAPACK.Scalar (zero)
import Numeric.LAPACK.Private (copyToColumnMajor, peekCInt, argMsg)

import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.Marshal.Alloc (alloca)
import Foreign.C.Types (CInt)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Ptr (Ptr)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)

import Text.Printf (printf)


solver ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width, Eq height,
    Class.Floating a) =>
   String -> height ->
   (Int -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr CInt -> ContT () IO ()) ->
   Full vert horiz height width a ->
   Full vert horiz height width a
solver name sh f (Array (MatrixShape.Full order extent) b) =
   Array.unsafeCreate (MatrixShape.Full ColumnMajor extent) $
      \xPtr -> do
   let (height,width) = Extent.dimensions extent
   Call.assert (name ++ ": height shapes mismatch") (sh == height)
   let n = Shape.size height
   let nrhs = Shape.size width
   evalContT $ do
      nPtr <- Call.cint n
      nrhsPtr <- Call.cint nrhs
      bPtr <- ContT $ withForeignPtr b
      liftIO $ copyToColumnMajor order n nrhs bPtr xPtr
      ldxPtr <- Call.leadingDim n
      f n nPtr nrhsPtr xPtr ldxPtr


withDeterminantInfo ::
   (Class.Floating a) =>
   String -> (Ptr CInt -> IO ()) -> IO a -> IO a
withDeterminantInfo name computation evaluation = alloca $ \infoPtr -> do
   computation infoPtr
   info <- peekCInt infoPtr
   case compare info (0::Int) of
      LT -> error $ printf argMsg name (-info)
      GT -> return zero
      EQ -> evaluation


withInfo :: String -> (Ptr CInt -> IO ()) -> IO ()
withInfo = Private.withInfo diagonalMsg

diagonalMsg :: String
diagonalMsg = "%d-th diagonal value is zero"
