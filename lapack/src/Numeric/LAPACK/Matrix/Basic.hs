{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Numeric.LAPACK.Matrix.Basic where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Matrix.Shape.Private (Order(RowMajor, ColumnMajor))
import Numeric.LAPACK.Matrix.Private (Full, General)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (zero, one)
import Numeric.LAPACK.Private (pointerSeq)

import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Storable (poke, peek)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)


transpose ::
   (Extent.C vert, Extent.C horiz) =>
   Full vert horiz height width a -> Full horiz vert width height a
transpose = Array.mapShape MatrixShape.transpose


singleRow :: Order -> Vector width a -> General () width a
singleRow order = Array.mapShape (MatrixShape.general order ())

singleColumn :: Order -> Vector height a -> General height () a
singleColumn order = Array.mapShape (flip (MatrixShape.general order) ())

flattenRow :: General () width a -> Vector width a
flattenRow = Array.mapShape MatrixShape.fullWidth

flattenColumn :: General height () a -> Vector height a
flattenColumn = Array.mapShape MatrixShape.fullHeight


forceRowMajor ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Full vert horiz height width a ->
   Full vert horiz height width a
forceRowMajor (Array shape@(MatrixShape.Full order extent) x) =
   case order of
      RowMajor -> Array shape x
      ColumnMajor ->
         Array.unsafeCreate (MatrixShape.Full RowMajor extent) $ \yPtr ->
         withForeignPtr x $ \xPtr -> do
            let (height, width) = Extent.dimensions extent
            let n = Shape.size width
            let m = Shape.size height
            Private.copyTransposed n m xPtr n yPtr

forceOrder ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Order ->
   Full vert horiz height width a ->
   Full vert horiz height width a
forceOrder order =
   case order of
      RowMajor -> forceRowMajor
      ColumnMajor -> transpose . forceRowMajor . transpose


scaleRows ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Class.Floating a) =>
   Vector height a ->
   Full vert horiz height width a ->
   Full vert horiz height width a
scaleRows
   (Array heightX x) (Array shape@(MatrixShape.Full order extent) a) =
      Array.unsafeCreate shape $ \bPtr -> do
   let (height,width) = Extent.dimensions extent
   Call.assert "scaleRows: sizes mismatch" (heightX == height)
   case order of
      RowMajor -> evalContT $ do
         let m = Shape.size height
         let n = Shape.size width
         alphaPtr <- Call.alloca
         nPtr <- Call.cint n
         xPtr <- ContT $ withForeignPtr x
         aPtr <- ContT $ withForeignPtr a
         incaPtr <- Call.cint 1
         incbPtr <- Call.cint 1
         liftIO $ sequence_ $ take m $
            zipWith3
               (\xkPtr akPtr bkPtr -> do
                  poke alphaPtr =<< peek xkPtr
                  BlasGen.copy nPtr akPtr incaPtr bkPtr incbPtr
                  BlasGen.scal nPtr alphaPtr bkPtr incbPtr)
               (pointerSeq 1 xPtr)
               (pointerSeq n aPtr)
               (pointerSeq n bPtr)
      ColumnMajor -> evalContT $ do
         let m = Shape.size width
         let n = Shape.size height
         transPtr <- Call.char 'N'
         nPtr <- Call.cint n
         klPtr <- Call.cint 0
         kuPtr <- Call.cint 0
         alphaPtr <- Call.number one
         xPtr <- ContT $ withForeignPtr x
         ldxPtr <- Call.leadingDim 1
         aPtr <- ContT $ withForeignPtr a
         incaPtr <- Call.cint 1
         betaPtr <- Call.number zero
         incbPtr <- Call.cint 1
         liftIO $ sequence_ $ take m $
            zipWith
               (\akPtr bkPtr ->
                  Private.gbmv transPtr
                     nPtr nPtr klPtr kuPtr alphaPtr xPtr ldxPtr
                     akPtr incaPtr betaPtr bkPtr incbPtr)
               (pointerSeq n aPtr)
               (pointerSeq n bPtr)

scaleColumns ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Eq width, Class.Floating a) =>
   Vector width a ->
   Full vert horiz height width a ->
   Full vert horiz height width a
scaleColumns x = transpose . scaleRows x . transpose
