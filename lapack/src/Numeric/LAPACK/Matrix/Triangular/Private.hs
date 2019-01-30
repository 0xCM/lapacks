{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.Triangular.Private where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor,ColumnMajor), flipOrder, uploFromOrder,
          Empty, Filled, NonUnit)
import Numeric.LAPACK.Matrix.Private (Full, Conjugation(Conjugated))
import Numeric.LAPACK.Scalar (zero)
import Numeric.LAPACK.Private
         (pointerSeq, copyBlock, copyCondConjugateToTemp,
          pokeCInt, fill, withInfo, errorCodeMsg)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (advancePtr)
import Foreign.C.Types (CInt)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Ptr (Ptr)
import Foreign.Storable (Storable)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)

import Data.Foldable (forM_)


diagonalPointers :: (Storable a) => Order -> Int -> Ptr a -> [Ptr a]
diagonalPointers order n aPtr =
   take n $ scanl advancePtr aPtr $
   case order of
      RowMajor -> iterate pred n
      ColumnMajor -> iterate succ 2

diagonalPointerPairs ::
   (Storable a, Storable b) =>
   Order -> Int -> Ptr a -> Ptr b -> [(Ptr a, Ptr b)]
diagonalPointerPairs order n aPtr bPtr =
   zip (pointerSeq 1 aPtr) $ diagonalPointers order n bPtr


columnMajorPointers ::
   (Storable a) => Int -> Ptr a -> Ptr a -> [(Int, ((Ptr a, Ptr a), Ptr a))]
columnMajorPointers n fullPtr packedPtr =
   let ds = iterate succ 1
   in  take n $ zip ds $
       zip
         (zip (pointerSeq 1 fullPtr) (pointerSeq n fullPtr))
         (scanl advancePtr packedPtr ds)

rowMajorPointers ::
   (Storable a) => Int -> Ptr a -> Ptr a -> [(Int, (Ptr a, Ptr a))]
rowMajorPointers n fullPtr packedPtr =
   let ds = iterate pred n
   in  take n $ zip ds $
       zip (pointerSeq (n+1) fullPtr) (scanl advancePtr packedPtr ds)


forPointers :: [(Int, a)] -> (Ptr CInt -> a -> IO ()) -> IO ()
forPointers xs act =
   alloca $ \nPtr ->
   forM_ xs $ \(d,ptrs) -> do
      pokeCInt nPtr d
      act nPtr ptrs


copyTriangleToTemp ::
   Class.Floating a =>
   Conjugation -> Order -> Int -> ForeignPtr a -> ContT r IO (Ptr a)
copyTriangleToTemp conj order =
   copyCondConjugateToTemp (order==RowMajor && conj==Conjugated)


unpackToTemp ::
   Storable a =>
   (Int -> Ptr a -> Ptr a -> IO ()) ->
   Int -> ForeignPtr a -> ContT r IO (Ptr a)
unpackToTemp f n a = do
   apPtr <- ContT $ withForeignPtr a
   aPtr <- Call.allocaArray (n*n)
   liftIO $ f n apPtr aPtr
   return aPtr


unpack :: Class.Floating a => Order -> Int -> Ptr a -> Ptr a -> IO ()
unpack order n packedPtr fullPtr =
   evalContT $ do
      uploPtr <- Call.char $ uploFromOrder order
      nPtr <- Call.cint n
      ldaPtr <- Call.leadingDim n
      liftIO $ withInfo errorCodeMsg "tpttr" $
         LapackGen.tpttr uploPtr nPtr packedPtr fullPtr ldaPtr

pack :: Class.Floating a => Order -> Int -> Ptr a -> Ptr a -> IO ()
pack order n = packRect order n n

packRect :: Class.Floating a => Order -> Int -> Int -> Ptr a -> Ptr a -> IO ()
packRect order n ld fullPtr packedPtr =
   evalContT $ do
      uploPtr <- Call.char $ uploFromOrder order
      nPtr <- Call.cint n
      ldaPtr <- Call.leadingDim ld
      liftIO $ withInfo errorCodeMsg "trttp" $
         LapackGen.trttp uploPtr nPtr fullPtr ldaPtr packedPtr


unpackZero, _unpackZero ::
   Class.Floating a => Order -> Int -> Ptr a -> Ptr a -> IO ()
_unpackZero order n packedPtr fullPtr = do
   fill zero (n*n) fullPtr
   unpack order n packedPtr fullPtr

unpackZero order n packedPtr fullPtr = do
   fillTriangle zero (flipOrder order) n fullPtr
   unpack order n packedPtr fullPtr

fillTriangle :: Class.Floating a => a -> Order -> Int -> Ptr a -> IO ()
fillTriangle z order n aPtr = evalContT $ do
   uploPtr <- Call.char $ uploFromOrder order
   nPtr <- Call.cint n
   zPtr <- Call.number z
   liftIO $ LapackGen.laset uploPtr nPtr nPtr zPtr zPtr aPtr nPtr


type Triangular lo diag up sh = Array (MatrixShape.Triangular lo diag up sh)

type FlexDiagonal diag sh =
         Triangular MatrixShape.Empty diag MatrixShape.Empty sh

newtype MultiplyRight diag sh a b lo up =
   MultiplyRight {getMultiplyRight :: Triangular lo diag up sh a -> b}

newtype Map diag sh a lo up =
   Map {getMap :: Triangular lo diag up sh a -> Triangular lo diag up sh a}

newtype Power diag sh a lo up =
   Power {
      getPower ::
         Triangular lo diag up sh a ->
         Triangular lo (PowerDiag lo up diag) up sh a
   }

type family PowerDiag lo up diag
type instance PowerDiag Empty up diag = diag
type instance PowerDiag Filled Empty diag = diag
type instance PowerDiag Filled Filled diag = NonUnit

caseTriDiagArray ::
   (MatrixShape.TriDiag diag) =>
   (Triangular lo diag up sh a -> b) ->
   (Triangular lo diag up sh a -> b) ->
   (Triangular lo diag up sh a -> b)
caseTriDiagArray fu fn a =
   MatrixShape.caseTriDiag
      (MatrixShape.triangularDiag $ Array.shape a)
      (fu a) (fn a)

multiplyDiagonal ::
   (Eq sh, MatrixShape.TriDiag diag) =>
   String ->
   (b -> sh) ->
   (Triangular lo diag up sh a -> b -> b) ->
   (Triangular lo diag up sh a -> b -> b)
multiplyDiagonal msg shape =
   caseTriDiagArray
      (\a b ->
         if MatrixShape.triangularSize (Array.shape a) == shape b
           then b
           else error ("Triangular." ++ msg))


fromBanded ::
   (Class.Floating a) =>
   Int -> Order -> Int -> ForeignPtr a -> Int -> Ptr a -> IO ()
fromBanded k order n a bSize bPtr =
   withForeignPtr a $ \aPtr -> do
      fill zero bSize bPtr
      let lda = k+1
      let pointers =
            zip [0..] $ zip (pointerSeq lda aPtr) $
            diagonalPointers order n bPtr
      case order of
         ColumnMajor ->
            forM_ pointers $ \(i,(xPtr,yPtr)) ->
               let j = min i k
               in copyBlock (j+1) (advancePtr xPtr (k-j)) (advancePtr yPtr (-j))
         RowMajor ->
            forM_ pointers $ \(i,(xPtr,yPtr)) ->
               copyBlock (min lda (n-i)) xPtr yPtr


type FlexLower diag sh = Array (MatrixShape.LowerTriangular diag sh)

takeLower ::
   (MatrixShape.TriDiag diag,
    Extent.C horiz, Shape.C height, Shape.C width, Class.Floating a) =>
   (diag, Order -> Int -> Ptr a -> IO ()) ->
   Full Extent.Small horiz height width a -> FlexLower diag height a
takeLower (diag, fillDiag) (Array (MatrixShape.Full order extent) a) =
   let (height,width) = Extent.dimensions extent
       m = Shape.size height
       n = Shape.size width
       k = case order of RowMajor -> n; ColumnMajor -> m
   in Array.unsafeCreate
         (MatrixShape.Triangular diag MatrixShape.lower order height) $ \lPtr ->
      withForeignPtr a $ \aPtr -> do
         let dstOrder = flipOrder order
         packRect dstOrder m k aPtr lPtr
         fillDiag dstOrder m lPtr
