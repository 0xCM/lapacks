{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Private where

import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor, ColumnMajor), transposeFromOrder)
import Numeric.LAPACK.Wrapper (Flip(Flip, getFlip))

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.LAPACK.FFI.Complex as LapackComplex
import qualified Numeric.BLAS.FFI.Real as BlasReal
import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class
import Numeric.LAPACK.Scalar (zero, one, isZero)

import qualified Foreign.Marshal.Utils as Marshal
import qualified Foreign.C.String as CStr
import Foreign.Marshal.Array (copyArray, advancePtr)
import Foreign.Marshal.Alloc (alloca)
import Foreign.C.Types (CChar, CInt)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Ptr (Ptr)
import Foreign.Storable (Storable, poke, peek)

import Text.Printf (printf)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT, runContT)
import Control.Monad.IO.Class (liftIO)
import Control.Monad (when, foldM)
import Control.Applicative ((<$>))

import qualified Data.Array.Comfort.Storable.Unchecked.Monadic as ArrayIO
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)

import qualified Data.Complex as Complex
import Data.Complex (Complex)
import Data.Tuple.HT (swap)

import Prelude hiding (sum)


fill :: (Class.Floating a) => a -> Int -> Ptr a -> IO ()
fill a n dstPtr = evalContT $ do
   nPtr <- Call.cint n
   srcPtr <- Call.number a
   incxPtr <- Call.cint 0
   incyPtr <- Call.cint 1
   liftIO $ BlasGen.copy nPtr srcPtr incxPtr dstPtr incyPtr


copyBlock :: (Class.Floating a) => Int -> Ptr a -> Ptr a -> IO ()
copyBlock n srcPtr dstPtr = evalContT $ do
   nPtr <- Call.cint n
   incxPtr <- Call.cint 1
   incyPtr <- Call.cint 1
   liftIO $ BlasGen.copy nPtr srcPtr incxPtr dstPtr incyPtr

copyToTemp :: (Storable a) => Int -> ForeignPtr a -> ContT r IO (Ptr a)
copyToTemp n fptr = do
   ptr <- ContT $ withForeignPtr fptr
   tmpPtr <- Call.allocaArray n
   liftIO $ copyArray tmpPtr ptr n
   return tmpPtr


{- |
Make a temporary copy only for complex matrices.
-}
conjugateToTemp ::
   (Class.Floating a) => Int -> ForeignPtr a -> ContT r IO (Ptr a)
conjugateToTemp n =
   runCopyToTemp $
   Class.switchFloating
      (CopyToTemp $ ContT . withForeignPtr)
      (CopyToTemp $ ContT . withForeignPtr)
      (CopyToTemp $ complexConjugateToTemp n)
      (CopyToTemp $ complexConjugateToTemp n)

newtype CopyToTemp r a =
   CopyToTemp {runCopyToTemp :: ForeignPtr a -> ContT r IO (Ptr a)}

complexConjugateToTemp ::
   Class.Real a =>
   Int -> ForeignPtr (Complex a) -> ContT r IO (Ptr (Complex a))
complexConjugateToTemp n x = do
   nPtr <- Call.cint n
   xPtr <- copyToTemp n x
   incxPtr <- Call.cint 1
   liftIO $ LapackComplex.lacgv nPtr xPtr incxPtr
   return xPtr


copyConjugate ::
   (Class.Floating a) =>
   Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
copyConjugate nPtr xPtr incxPtr yPtr incyPtr = do
   BlasGen.copy nPtr xPtr incxPtr yPtr incyPtr
   lacgv nPtr yPtr incyPtr

copyCondConjugate ::
   (Class.Floating a) =>
   Bool -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
copyCondConjugate conj nPtr xPtr incxPtr yPtr incyPtr = do
   BlasGen.copy nPtr xPtr incxPtr yPtr incyPtr
   when conj $ lacgv nPtr yPtr incyPtr

condConjugateToTemp ::
   (Class.Floating a) =>
   Bool -> Int -> ForeignPtr a -> ContT r IO (Ptr a)
condConjugateToTemp conj n x =
   if conj then conjugateToTemp n x else ContT $ withForeignPtr x

copyCondConjugateToTemp ::
   (Class.Floating a) =>
   Bool -> Int -> ForeignPtr a -> ContT r IO (Ptr a)
copyCondConjugateToTemp conj n a = do
   bPtr <- Call.allocaArray n
   liftIO $ evalContT $ do
      aPtr <- ContT $ withForeignPtr a
      sizePtr <- Call.cint n
      incPtr <- Call.cint 1
      liftIO $ copyCondConjugate conj sizePtr aPtr incPtr bPtr incPtr
      return bPtr



{- |
In ColumnMajor:
Copy a m-by-n-matrix with lda>=m and ldb>=m.
-}
copySubMatrix ::
   (Class.Floating a) =>
   Int -> Int -> Int -> Ptr a -> Int -> Ptr a -> IO ()
copySubMatrix = copySubTrapezoid 'A'

copySubTrapezoid ::
   (Class.Floating a) =>
   Char -> Int -> Int -> Int -> Ptr a -> Int -> Ptr a -> IO ()
copySubTrapezoid side m n lda aPtr ldb bPtr = evalContT $ do
   uploPtr <- Call.char side
   mPtr <- Call.cint m
   nPtr <- Call.cint n
   ldaPtr <- Call.leadingDim lda
   ldbPtr <- Call.leadingDim ldb
   liftIO $ LapackGen.lacpy uploPtr mPtr nPtr aPtr ldaPtr bPtr ldbPtr

copyTransposed ::
   (Class.Floating a) =>
   Int -> Int -> Ptr a -> Int -> Ptr a -> IO ()
copyTransposed n m aPtr ldb bPtr = evalContT $ do
   nPtr <- Call.cint n
   incaPtr <- Call.cint m
   incbPtr <- Call.cint 1
   liftIO $ sequence_ $ take m $
      zipWith
         (\akPtr bkPtr -> BlasGen.copy nPtr akPtr incaPtr bkPtr incbPtr)
         (pointerSeq 1 aPtr)
         (pointerSeq ldb bPtr)


{- |
Copy a m-by-n-matrix to ColumnMajor order.
-}
copyToColumnMajor ::
   (Class.Floating a) =>
   Order -> Int -> Int -> Ptr a -> Ptr a -> IO ()
copyToColumnMajor order m n aPtr bPtr =
   case order of
      RowMajor -> copyTransposed m n aPtr m bPtr
      ColumnMajor -> copyBlock (m*n) aPtr bPtr

copyToSubColumnMajor ::
   (Class.Floating a) =>
   Order -> Int -> Int -> Ptr a -> Int -> Ptr a -> IO ()
copyToSubColumnMajor order m n aPtr ldb bPtr =
   case order of
      RowMajor -> copyTransposed m n aPtr ldb bPtr
      ColumnMajor ->
         if m==ldb
           then copyBlock (m*n) aPtr bPtr
           else copySubMatrix m n m aPtr ldb bPtr


pointerSeq :: (Storable a) => Int -> Ptr a -> [Ptr a]
pointerSeq k ptr = iterate (flip advancePtr k) ptr


createHigherArray ::
   (Shape.C sh, Class.Floating a) =>
   sh -> Int -> Int -> Int ->
   ((Ptr a, Int) -> IO rank) -> IO (rank, Array sh a)
createHigherArray shapeX m n nrhs act =
   fmap swap $ ArrayIO.unsafeCreateWithSizeAndResult shapeX $ \ _ xPtr ->
   if m>n
      then
         runContT (Call.allocaArray (m*nrhs)) $ \tmpPtr -> do
            r <- act (tmpPtr,m)
            copySubMatrix n nrhs m tmpPtr n xPtr
            return r
      else act (xPtr,n)



newtype Sum a = Sum {runSum :: Int -> Ptr a -> Int -> IO a}

sum :: Class.Floating a => Int -> Ptr a -> Int -> IO a
sum =
   runSum $
   Class.switchFloating
      (Sum sumReal)
      (Sum sumReal)
      (Sum sumComplex)
      (Sum sumComplex)

sumReal :: Class.Real a => Int -> Ptr a -> Int -> IO a
sumReal n xPtr incx =
   evalContT $ do
      nPtr <- Call.cint n
      incxPtr <- Call.cint incx
      yPtr <- Call.real one
      incyPtr <- Call.cint 0
      liftIO $ BlasReal.dot nPtr xPtr incxPtr yPtr incyPtr

sumComplex :: Class.Real a => Int -> Ptr (Complex a) -> Int -> IO (Complex a)
sumComplex n xPtr incx =
   evalContT $ do
      transPtr <- Call.char 'N'
      mPtr <- Call.cint 1
      nPtr <- Call.cint n
      alphaPtr <- Call.number one
      onePtr <- Call.number one
      zeroincPtr <- Call.cint 0
      aPtr <- Call.allocaArray n
      ldaPtr <- Call.leadingDim 1
      incxPtr <- Call.cint incx
      betaPtr <- Call.number zero
      yPtr <- Call.alloca
      incyPtr <- Call.cint 1
      liftIO $ do
         BlasGen.copy nPtr onePtr zeroincPtr aPtr incyPtr
         gemv
            transPtr mPtr nPtr alphaPtr aPtr ldaPtr
            xPtr incxPtr betaPtr yPtr incyPtr
         peek yPtr


product :: Class.Floating a => Int -> Ptr a -> Int -> IO a
product n xPtr incx =
   foldM (\x ptr -> do y <- peek ptr; return $! x*y) one $
   take n $ pointerSeq incx xPtr


newtype LACGV a = LACGV {getLACGV :: Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

lacgv :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
lacgv =
   getLACGV $
   Class.switchFloating
      (LACGV $ const $ const $ const $ return ())
      (LACGV $ const $ const $ const $ return ())
      (LACGV LapackComplex.lacgv)
      (LACGV LapackComplex.lacgv)


{-
Work around an inconsistency of BLAS.
In case of a zero-column matrix
BLAS's gemv and gbmv do not initialize the target vector.
In contrast, these work-arounds do.
-}
{-# INLINE gemv #-}
gemv ::
   (Class.Floating a) =>
   Ptr CChar -> Ptr CInt -> Ptr CInt ->
   Ptr a -> Ptr a -> Ptr CInt ->
   Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
gemv transPtr mPtr nPtr
      alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr = do
   initializeMV transPtr mPtr nPtr betaPtr yPtr incyPtr
   BlasGen.gemv transPtr mPtr nPtr
      alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

{-# INLINE gbmv #-}
gbmv ::
   (Class.Floating a) =>
   Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CInt ->
   Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt ->
   Ptr a -> Ptr a -> Ptr CInt -> IO ()
gbmv transPtr mPtr nPtr klPtr kuPtr
      alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr = do
   initializeMV transPtr mPtr nPtr betaPtr yPtr incyPtr
   BlasGen.gbmv transPtr mPtr nPtr klPtr kuPtr
      alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

initializeMV ::
   Class.Floating a =>
   Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
initializeMV transPtr mPtr nPtr betaPtr yPtr incyPtr = do
   trans <- peek transPtr
   let (mtPtr,ntPtr) =
         if trans == CStr.castCharToCChar 'N'
            then (mPtr,nPtr) else (nPtr,mPtr)
   n <- peek ntPtr
   beta <- peek betaPtr
   when (n == 0 && isZero beta) $
      Marshal.with 0 $ \incbPtr ->
      BlasGen.copy mtPtr betaPtr incbPtr yPtr incyPtr


multiplyMatrix ::
   (Class.Floating a) =>
   Order -> Order -> Int -> Int -> Int ->
   ForeignPtr a -> ForeignPtr a -> Ptr a -> IO ()
multiplyMatrix orderA orderB m k n a b cPtr = do
   let lda = case orderA of RowMajor -> k; ColumnMajor -> m
   let ldb = case orderB of RowMajor -> n; ColumnMajor -> k
   let ldc = m
   evalContT $ do
      transaPtr <- Call.char $ transposeFromOrder orderA
      transbPtr <- Call.char $ transposeFromOrder orderB
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim lda
      bPtr <- ContT $ withForeignPtr b
      ldbPtr <- Call.leadingDim ldb
      betaPtr <- Call.number zero
      ldcPtr <- Call.leadingDim ldc
      liftIO $
         BlasGen.gemm
            transaPtr transbPtr mPtr nPtr kPtr alphaPtr aPtr ldaPtr
            bPtr ldbPtr betaPtr cPtr ldcPtr



withAutoWorkspaceInfo ::
   (Class.Floating a) =>
   String -> String -> (Ptr a -> Ptr CInt -> Ptr CInt -> IO ()) -> IO ()
withAutoWorkspaceInfo msg name computation =
   withInfo msg name $ \infoPtr ->
   withAutoWorkspace $ \workPtr lworkPtr ->
      computation workPtr lworkPtr infoPtr

withAutoWorkspace ::
   (Class.Floating a) =>
   (Ptr a -> Ptr CInt -> IO ()) -> IO ()
withAutoWorkspace computation = evalContT $ do
   lworkPtr <- Call.cint (-1)
   lwork <- liftIO $ alloca $ \workPtr -> do
      computation workPtr lworkPtr
      max 1 . ceilingSize <$> peek workPtr
   workPtr <- Call.allocaArray lwork
   liftIO $ pokeCInt lworkPtr lwork
   liftIO $ computation workPtr lworkPtr

withInfo :: String -> String -> (Ptr CInt -> IO ()) -> IO ()
withInfo msg name computation = alloca $ \infoPtr -> do
   computation infoPtr
   info <- peekCInt infoPtr
   case compare info (0::Int) of
      EQ -> return ()
      LT -> error $ printf argMsg name (-info)
      GT -> error $ name ++ ": " ++ printf msg info

argMsg :: String
argMsg = "%s: illegal value in %d-th argument"

errorCodeMsg :: String
errorCodeMsg = "unknown error code %d"

rankMsg :: String
rankMsg = "deficient rank %d"

definiteMsg :: String
definiteMsg = "minor of order %d not positive definite"

eigenMsg :: String
eigenMsg = "%d off-diagonal elements not converging"


pokeCInt :: Ptr CInt -> Int -> IO ()
pokeCInt ptr = poke ptr . fromIntegral

peekCInt :: Ptr CInt -> IO Int
peekCInt ptr = fromIntegral <$> peek ptr


ceilingSize :: (Class.Floating a) => a -> Int
ceilingSize =
   getFlip $
   Class.switchFloating
      (Flip ceiling)
      (Flip ceiling)
      (Flip $ ceiling . Complex.realPart)
      (Flip $ ceiling . Complex.realPart)
