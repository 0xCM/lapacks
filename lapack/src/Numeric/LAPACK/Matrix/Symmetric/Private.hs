module Numeric.LAPACK.Matrix.Symmetric.Private where

import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Triangular.Private
         (diagonalPointerPairs, columnMajorPointers, rowMajorPointers,
          forPointers, pack, unpackToTemp, copyTriangleToTemp)
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor,ColumnMajor), uploFromOrder)
import Numeric.LAPACK.Matrix.Private
         (Full, Conjugation(NonConjugated, Conjugated))
import Numeric.LAPACK.Linear.Private (solver, withDeterminantInfo, withInfo)
import Numeric.LAPACK.Scalar (zero, one)
import Numeric.LAPACK.Private (copyBlock, copyToTemp, copyCondConjugate)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Shape (triangleSize)

import Foreign.Marshal.Array (advancePtr)
import Foreign.C.Types (CInt)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Ptr (Ptr)
import Foreign.Storable (Storable, peek)

import qualified System.IO.Lazy as LazyIO

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Applicative ((<$>))


unpack :: Class.Floating a =>
   Conjugation -> Order -> Int -> Ptr a -> Ptr a -> IO ()
unpack conj order n packedPtr fullPtr = evalContT $ do
   incxPtr <- Call.cint 1
   incyPtr <- Call.cint n
   liftIO $ case order of
      RowMajor ->
         forPointers (rowMajorPointers n fullPtr packedPtr) $
               \nPtr (dstPtr,srcPtr) -> do
            copyCondConjugate (conj==Conjugated)
               nPtr srcPtr incxPtr dstPtr incyPtr
            BlasGen.copy nPtr srcPtr incxPtr dstPtr incxPtr
      ColumnMajor ->
         forPointers (columnMajorPointers n fullPtr packedPtr) $
               \nPtr ((dstRowPtr,dstColumnPtr),srcPtr) -> do
            copyCondConjugate (conj==Conjugated)
               nPtr srcPtr incxPtr dstRowPtr incyPtr
            BlasGen.copy nPtr srcPtr incxPtr dstColumnPtr incxPtr


square ::
   (Class.Floating a) =>
   Conjugation -> Order -> Int -> ForeignPtr a -> Ptr a -> IO ()
square conj order n a bpPtr =
   evalContT $ do
      sidePtr <- Call.char 'L'
      uploPtr <- Call.char 'U'
      nPtr <- Call.cint n
      ldPtr <- Call.leadingDim n
      aPtr <- unpackToTemp (unpack conj order) n a
      bPtr <- Call.allocaArray (n*n)
      alphaPtr <- Call.number one
      betaPtr <- Call.number zero
      liftIO $ do
         (if conj==Conjugated then BlasGen.hemm else BlasGen.symm)
            sidePtr uploPtr
            nPtr nPtr alphaPtr aPtr ldPtr
            aPtr ldPtr betaPtr bPtr ldPtr
         pack order n bPtr bpPtr


solve ::
   (Extent.C vert, Extent.C horiz,
    Shape.C width, Shape.C height, Eq height, Class.Floating a) =>
   String -> Conjugation -> Order -> height -> ForeignPtr a ->
   Full vert horiz height width a ->
   Full vert horiz height width a
solve name conj order sh a =
   solver name sh $ \n nPtr nrhsPtr xPtr ldxPtr -> do
      uploPtr <- Call.char $ uploFromOrder order
      apPtr <- copyTriangleToTemp conj order (triangleSize n) a
      ipivPtr <- Call.allocaArray n
      liftIO $
         let (lapackName,slv) =
               case conj of
                  Conjugated -> ("hpsv", LapackGen.hpsv)
                  NonConjugated -> ("spsv", LapackGen.spsv)
         in withInfo lapackName $
               slv uploPtr nPtr nrhsPtr apPtr ipivPtr xPtr ldxPtr


inverse ::
   Class.Floating a =>
   Conjugation -> Order -> Int -> ForeignPtr a -> Int -> Ptr a -> IO ()
inverse conj order n a triSize bPtr = evalContT $ do
   uploPtr <- Call.char $ uploFromOrder order
   nPtr <- Call.cint n
   aPtr <- ContT $ withForeignPtr a
   ipivPtr <- Call.allocaArray n
   workPtr <- Call.allocaArray n
   liftIO $ do
      copyBlock triSize aPtr bPtr
      case conj of
         Conjugated -> do
            withInfo "hptrf" $ LapackGen.hptrf uploPtr nPtr bPtr ipivPtr
            withInfo "hptri" $ LapackGen.hptri uploPtr nPtr bPtr ipivPtr workPtr
         NonConjugated -> do
            withInfo "sptrf" $ LapackGen.sptrf uploPtr nPtr bPtr ipivPtr
            withInfo "sptri" $ LapackGen.sptri uploPtr nPtr bPtr ipivPtr workPtr


blockDiagonalPointers ::
   (Storable a) =>
   Order -> [(Ptr CInt, Ptr a)] -> LazyIO.T [(Ptr a, Maybe (Ptr a, Ptr a))]
blockDiagonalPointers order =
   let go ((ipiv0Ptr,a0Ptr):ptrs0) = do
         ipiv <- LazyIO.interleave $ peek ipiv0Ptr
         (ext,ptrTuples) <-
            if ipiv >= 0
               then (,) Nothing <$> go ptrs0
               else
                  case ptrs0 of
                     [] -> error "Symmetric.determinant: incomplete 2x2 block"
                     (_ipiv1Ptr,a1Ptr):ptrs1 ->
                        let bPtr =
                              case order of
                                 ColumnMajor -> advancePtr a1Ptr (-1)
                                 RowMajor -> advancePtr a0Ptr 1
                        in (,) (Just (a1Ptr,bPtr)) <$> go ptrs1
         return $ (a0Ptr,ext) : ptrTuples
       go [] = return []
   in go

determinant ::
   (Class.Floating a, Class.Floating ar) =>
   Conjugation -> ((Ptr a, Maybe (Ptr a, Ptr a)) -> IO ar) ->
   Order -> Int -> ForeignPtr a -> IO ar
determinant conj peekBlockDeterminant order n a = evalContT $ do
   uploPtr <- Call.char $ uploFromOrder order
   nPtr <- Call.cint n
   aPtr <- copyToTemp (triangleSize n) a
   ipivPtr <- Call.allocaArray n
   let (name,trf) =
         case conj of
            Conjugated -> ("hptrf", LapackGen.hptrf)
            NonConjugated -> ("sptrf", LapackGen.sptrf)
   liftIO $ withDeterminantInfo name
      (trf uploPtr nPtr aPtr ipivPtr)
      (((return $!) =<<) $
       LazyIO.run
         (fmap product $
          mapM (LazyIO.interleave . peekBlockDeterminant) =<<
          blockDiagonalPointers order
            (diagonalPointerPairs order n ipivPtr aPtr)))
