-- Do not edit! Automatically generated by create-lapack-ffi.
module Numeric.BLAS.CArray.ComplexDouble where

import qualified Numeric.BLAS.FFI.ComplexDouble as FFI
import qualified Numeric.Netlib.CArray.Utility as Call

import Data.Array.IOCArray (IOCArray, getBounds)
import Data.Array.CArray (CArray, bounds)

import Data.Complex (Complex)

import Foreign.Storable.Complex ()
import Foreign.Storable (peek)
import Foreign.C.Types (CInt)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Applicative (pure, (<$>), (<*>))


axpy ::
   Int {- ^ n -} ->
   Complex Double {- ^ za -} ->
   CArray Int (Complex Double) {- ^ zx -} ->
   Int {- ^ incx -} ->
   IOCArray Int (Complex Double) {- ^ zy -} ->
   Int {- ^ incy -} ->
   IO ()
axpy n za zx incx zy incy = do
   let zxDim0 = Call.sizes1 $ bounds zx
   zyDim0 <- Call.sizes1 <$> getBounds zy
   Call.assert "axpy: 1+(n-1)*abs(incx) == zxDim0" (1+(n-1)*abs(incx) == zxDim0)
   Call.assert "axpy: 1+(n-1)*abs(incy) == zyDim0" (1+(n-1)*abs(incy) == zyDim0)
   evalContT $ do
      nPtr <- Call.cint n
      zaPtr <- Call.complexDouble za
      zxPtr <- Call.array zx
      incxPtr <- Call.cint incx
      zyPtr <- Call.ioarray zy
      incyPtr <- Call.cint incy
      liftIO $ FFI.axpy nPtr zaPtr zxPtr incxPtr zyPtr incyPtr

casum ::
   Int {- ^ n -} ->
   IOCArray Int (Complex Double) {- ^ zx -} ->
   Int {- ^ incx -} ->
   IO Double
casum n zx incx = do
   zxDim0 <- Call.sizes1 <$> getBounds zx
   Call.assert "casum: 1+(n-1)*abs(incx) == zxDim0" (1+(n-1)*abs(incx) == zxDim0)
   evalContT $ do
      nPtr <- Call.cint n
      zxPtr <- Call.ioarray zx
      incxPtr <- Call.cint incx
      liftIO $ FFI.casum nPtr zxPtr incxPtr

cnrm2 ::
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   IO Double
cnrm2 x incx = do
   let xDim0 = Call.sizes1 $ bounds x
   let n = xDim0
   evalContT $ do
      nPtr <- Call.cint n
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      liftIO $ FFI.cnrm2 nPtr xPtr incxPtr

copy ::
   Int {- ^ n -} ->
   CArray Int (Complex Double) {- ^ zx -} ->
   Int {- ^ incx -} ->
   Int {- ^ incy -} ->
   IO (CArray Int (Complex Double))
copy n zx incx incy = do
   let zxDim0 = Call.sizes1 $ bounds zx
   Call.assert "copy: 1+(n-1)*abs(incx) == zxDim0" (1+(n-1)*abs(incx) == zxDim0)
   zy <- Call.newArray1 (1+(n-1)*abs(incy))
   evalContT $ do
      nPtr <- Call.cint n
      zxPtr <- Call.array zx
      incxPtr <- Call.cint incx
      zyPtr <- Call.ioarray zy
      incyPtr <- Call.cint incy
      liftIO $ FFI.copy nPtr zxPtr incxPtr zyPtr incyPtr
      liftIO $ Call.freezeArray zy

gbmv ::
   Char {- ^ trans -} ->
   Int {- ^ m -} ->
   Int {- ^ kl -} ->
   Int {- ^ ku -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Double {- ^ beta -} ->
   IOCArray Int (Complex Double) {- ^ y -} ->
   Int {- ^ incy -} ->
   IO ()
gbmv trans m kl ku alpha a x incx beta y incy = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   let xDim0 = Call.sizes1 $ bounds x
   yDim0 <- Call.sizes1 <$> getBounds y
   let n = aDim0
   let lda = aDim1
   let _xSize = xDim0
   let _ySize = yDim0
   evalContT $ do
      transPtr <- Call.char trans
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      klPtr <- Call.cint kl
      kuPtr <- Call.cint ku
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexDouble beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.gbmv transPtr mPtr nPtr klPtr kuPtr alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

gemm ::
   Char {- ^ transa -} ->
   Char {- ^ transb -} ->
   Int {- ^ m -} ->
   Int {- ^ k -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   CArray (Int,Int) (Complex Double) {- ^ b -} ->
   Complex Double {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ c -} ->
   IO ()
gemm transa transb m k alpha a b beta c = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   let (bDim0,bDim1) = Call.sizes2 $ bounds b
   (cDim0,cDim1) <- Call.sizes2 <$> getBounds c
   let _ka = aDim0
   let lda = aDim1
   let _kb = bDim0
   let ldb = bDim1
   let n = cDim0
   let ldc = cDim1
   evalContT $ do
      transaPtr <- Call.char transa
      transbPtr <- Call.char transb
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.complexDouble beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.gemm transaPtr transbPtr mPtr nPtr kPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

gemv ::
   Char {- ^ trans -} ->
   Int {- ^ m -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Double {- ^ beta -} ->
   IOCArray Int (Complex Double) {- ^ y -} ->
   Int {- ^ incy -} ->
   IO ()
gemv trans m alpha a x incx beta y incy = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   let xDim0 = Call.sizes1 $ bounds x
   yDim0 <- Call.sizes1 <$> getBounds y
   let n = aDim0
   let lda = aDim1
   let _xSize = xDim0
   let _ySize = yDim0
   evalContT $ do
      transPtr <- Call.char trans
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexDouble beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.gemv transPtr mPtr nPtr alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

gerc ::
   Int {- ^ m -} ->
   Complex Double {- ^ alpha -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   CArray Int (Complex Double) {- ^ y -} ->
   Int {- ^ incy -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ a -} ->
   IO ()
gerc m alpha x incx y incy a = do
   let xDim0 = Call.sizes1 $ bounds x
   let yDim0 = Call.sizes1 $ bounds y
   (aDim0,aDim1) <- Call.sizes2 <$> getBounds a
   let _xSize = xDim0
   let _ySize = yDim0
   let n = aDim0
   let lda = aDim1
   evalContT $ do
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      yPtr <- Call.array y
      incyPtr <- Call.cint incy
      aPtr <- Call.ioarray a
      ldaPtr <- Call.cint lda
      liftIO $ FFI.gerc mPtr nPtr alphaPtr xPtr incxPtr yPtr incyPtr aPtr ldaPtr

geru ::
   Int {- ^ m -} ->
   Complex Double {- ^ alpha -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   CArray Int (Complex Double) {- ^ y -} ->
   Int {- ^ incy -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ a -} ->
   IO ()
geru m alpha x incx y incy a = do
   let xDim0 = Call.sizes1 $ bounds x
   let yDim0 = Call.sizes1 $ bounds y
   (aDim0,aDim1) <- Call.sizes2 <$> getBounds a
   let _xSize = xDim0
   let _ySize = yDim0
   let n = aDim0
   let lda = aDim1
   evalContT $ do
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      yPtr <- Call.array y
      incyPtr <- Call.cint incy
      aPtr <- Call.ioarray a
      ldaPtr <- Call.cint lda
      liftIO $ FFI.geru mPtr nPtr alphaPtr xPtr incxPtr yPtr incyPtr aPtr ldaPtr

hbmv ::
   Char {- ^ uplo -} ->
   Int {- ^ k -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Double {- ^ beta -} ->
   IOCArray Int (Complex Double) {- ^ y -} ->
   Int {- ^ incy -} ->
   IO ()
hbmv uplo k alpha a x incx beta y incy = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   let xDim0 = Call.sizes1 $ bounds x
   yDim0 <- Call.sizes1 <$> getBounds y
   let n = aDim0
   let lda = aDim1
   let _xSize = xDim0
   let _ySize = yDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexDouble beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.hbmv uploPtr nPtr kPtr alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

hemm ::
   Char {- ^ side -} ->
   Char {- ^ uplo -} ->
   Int {- ^ m -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   CArray (Int,Int) (Complex Double) {- ^ b -} ->
   Complex Double {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ c -} ->
   IO ()
hemm side uplo m alpha a b beta c = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   let (bDim0,bDim1) = Call.sizes2 $ bounds b
   (cDim0,cDim1) <- Call.sizes2 <$> getBounds c
   let _ka = aDim0
   let lda = aDim1
   let n = bDim0
   let ldb = bDim1
   let ldc = cDim1
   Call.assert "hemm: n == cDim0" (n == cDim0)
   evalContT $ do
      sidePtr <- Call.char side
      uploPtr <- Call.char uplo
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.complexDouble beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.hemm sidePtr uploPtr mPtr nPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

hemv ::
   Char {- ^ uplo -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Double {- ^ beta -} ->
   IOCArray Int (Complex Double) {- ^ y -} ->
   Int {- ^ incy -} ->
   IO ()
hemv uplo alpha a x incx beta y incy = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   let xDim0 = Call.sizes1 $ bounds x
   yDim0 <- Call.sizes1 <$> getBounds y
   let n = aDim0
   let lda = aDim1
   let _xSize = xDim0
   let _ySize = yDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexDouble beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.hemv uploPtr nPtr alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

her ::
   Char {- ^ uplo -} ->
   Double {- ^ alpha -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ a -} ->
   IO ()
her uplo alpha x incx a = do
   let xDim0 = Call.sizes1 $ bounds x
   (aDim0,aDim1) <- Call.sizes2 <$> getBounds a
   let _xSize = xDim0
   let n = aDim0
   let lda = aDim1
   evalContT $ do
      uploPtr <- Call.char uplo
      nPtr <- Call.cint n
      alphaPtr <- Call.double alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      aPtr <- Call.ioarray a
      ldaPtr <- Call.cint lda
      liftIO $ FFI.her uploPtr nPtr alphaPtr xPtr incxPtr aPtr ldaPtr

her2 ::
   Char {- ^ uplo -} ->
   Complex Double {- ^ alpha -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   CArray Int (Complex Double) {- ^ y -} ->
   Int {- ^ incy -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ a -} ->
   IO ()
her2 uplo alpha x incx y incy a = do
   let xDim0 = Call.sizes1 $ bounds x
   let yDim0 = Call.sizes1 $ bounds y
   (aDim0,aDim1) <- Call.sizes2 <$> getBounds a
   let _xSize = xDim0
   let _ySize = yDim0
   let n = aDim0
   let lda = aDim1
   evalContT $ do
      uploPtr <- Call.char uplo
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      yPtr <- Call.array y
      incyPtr <- Call.cint incy
      aPtr <- Call.ioarray a
      ldaPtr <- Call.cint lda
      liftIO $ FFI.her2 uploPtr nPtr alphaPtr xPtr incxPtr yPtr incyPtr aPtr ldaPtr

her2k ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Int {- ^ k -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   CArray (Int,Int) (Complex Double) {- ^ b -} ->
   Double {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ c -} ->
   IO ()
her2k uplo trans k alpha a b beta c = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   let (bDim0,bDim1) = Call.sizes2 $ bounds b
   (cDim0,cDim1) <- Call.sizes2 <$> getBounds c
   let _ka = aDim0
   let lda = aDim1
   let _kb = bDim0
   let ldb = bDim1
   let n = cDim0
   let ldc = cDim1
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.double beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.her2k uploPtr transPtr nPtr kPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

herk ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Int {- ^ k -} ->
   Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   Double {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ c -} ->
   IO ()
herk uplo trans k alpha a beta c = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   (cDim0,cDim1) <- Call.sizes2 <$> getBounds c
   let _ka = aDim0
   let lda = aDim1
   let n = cDim0
   let ldc = cDim1
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.double alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      betaPtr <- Call.double beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.herk uploPtr transPtr nPtr kPtr alphaPtr aPtr ldaPtr betaPtr cPtr ldcPtr

hpmv ::
   Char {- ^ uplo -} ->
   Int {- ^ n -} ->
   Complex Double {- ^ alpha -} ->
   CArray Int (Complex Double) {- ^ ap -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Double {- ^ beta -} ->
   IOCArray Int (Complex Double) {- ^ y -} ->
   Int {- ^ incy -} ->
   IO ()
hpmv uplo n alpha ap x incx beta y incy = do
   let apDim0 = Call.sizes1 $ bounds ap
   let xDim0 = Call.sizes1 $ bounds x
   yDim0 <- Call.sizes1 <$> getBounds y
   let _apSize = apDim0
   let _xSize = xDim0
   let _ySize = yDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      apPtr <- Call.array ap
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexDouble beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.hpmv uploPtr nPtr alphaPtr apPtr xPtr incxPtr betaPtr yPtr incyPtr

hpr ::
   Char {- ^ uplo -} ->
   Int {- ^ n -} ->
   Double {- ^ alpha -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   IOCArray Int (Complex Double) {- ^ ap -} ->
   IO ()
hpr uplo n alpha x incx ap = do
   let xDim0 = Call.sizes1 $ bounds x
   apDim0 <- Call.sizes1 <$> getBounds ap
   let _xSize = xDim0
   let _apSize = apDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      nPtr <- Call.cint n
      alphaPtr <- Call.double alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      apPtr <- Call.ioarray ap
      liftIO $ FFI.hpr uploPtr nPtr alphaPtr xPtr incxPtr apPtr

hpr2 ::
   Char {- ^ uplo -} ->
   Int {- ^ n -} ->
   Complex Double {- ^ alpha -} ->
   CArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   CArray Int (Complex Double) {- ^ y -} ->
   Int {- ^ incy -} ->
   IOCArray Int (Complex Double) {- ^ ap -} ->
   IO ()
hpr2 uplo n alpha x incx y incy ap = do
   let xDim0 = Call.sizes1 $ bounds x
   let yDim0 = Call.sizes1 $ bounds y
   apDim0 <- Call.sizes1 <$> getBounds ap
   let _xSize = xDim0
   let _ySize = yDim0
   let _apSize = apDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      yPtr <- Call.array y
      incyPtr <- Call.cint incy
      apPtr <- Call.ioarray ap
      liftIO $ FFI.hpr2 uploPtr nPtr alphaPtr xPtr incxPtr yPtr incyPtr apPtr

iamax ::
   Int {- ^ n -} ->
   CArray Int (Complex Double) {- ^ zx -} ->
   Int {- ^ incx -} ->
   IO CInt
iamax n zx incx = do
   let zxDim0 = Call.sizes1 $ bounds zx
   Call.assert "iamax: 1+(n-1)*abs(incx) == zxDim0" (1+(n-1)*abs(incx) == zxDim0)
   evalContT $ do
      nPtr <- Call.cint n
      zxPtr <- Call.array zx
      incxPtr <- Call.cint incx
      liftIO $ FFI.iamax nPtr zxPtr incxPtr

rotg ::
   Complex Double {- ^ ca -} ->
   Complex Double {- ^ cb -} ->
   IO (Double, Complex Double)
rotg ca cb = do
   evalContT $ do
      caPtr <- Call.complexDouble ca
      cbPtr <- Call.complexDouble cb
      cPtr <- Call.alloca
      sPtr <- Call.alloca
      liftIO $ FFI.rotg caPtr cbPtr cPtr sPtr
      liftIO $ pure (,)
         <*> peek cPtr
         <*> peek sPtr

rrot ::
   Int {- ^ n -} ->
   IOCArray Int (Complex Double) {- ^ cx -} ->
   Int {- ^ incx -} ->
   IOCArray Int (Complex Double) {- ^ cy -} ->
   Int {- ^ incy -} ->
   Double {- ^ c -} ->
   Double {- ^ s -} ->
   IO ()
rrot n cx incx cy incy c s = do
   cxDim0 <- Call.sizes1 <$> getBounds cx
   cyDim0 <- Call.sizes1 <$> getBounds cy
   let _cxSize = cxDim0
   let _cySize = cyDim0
   evalContT $ do
      nPtr <- Call.cint n
      cxPtr <- Call.ioarray cx
      incxPtr <- Call.cint incx
      cyPtr <- Call.ioarray cy
      incyPtr <- Call.cint incy
      cPtr <- Call.double c
      sPtr <- Call.double s
      liftIO $ FFI.rrot nPtr cxPtr incxPtr cyPtr incyPtr cPtr sPtr

rscal ::
   Int {- ^ n -} ->
   Double {- ^ da -} ->
   IOCArray Int (Complex Double) {- ^ zx -} ->
   Int {- ^ incx -} ->
   IO ()
rscal n da zx incx = do
   zxDim0 <- Call.sizes1 <$> getBounds zx
   Call.assert "rscal: 1+(n-1)*abs(incx) == zxDim0" (1+(n-1)*abs(incx) == zxDim0)
   evalContT $ do
      nPtr <- Call.cint n
      daPtr <- Call.double da
      zxPtr <- Call.ioarray zx
      incxPtr <- Call.cint incx
      liftIO $ FFI.rscal nPtr daPtr zxPtr incxPtr

scal ::
   Int {- ^ n -} ->
   Complex Double {- ^ za -} ->
   IOCArray Int (Complex Double) {- ^ zx -} ->
   Int {- ^ incx -} ->
   IO ()
scal n za zx incx = do
   zxDim0 <- Call.sizes1 <$> getBounds zx
   Call.assert "scal: 1+(n-1)*abs(incx) == zxDim0" (1+(n-1)*abs(incx) == zxDim0)
   evalContT $ do
      nPtr <- Call.cint n
      zaPtr <- Call.complexDouble za
      zxPtr <- Call.ioarray zx
      incxPtr <- Call.cint incx
      liftIO $ FFI.scal nPtr zaPtr zxPtr incxPtr

swap ::
   Int {- ^ n -} ->
   IOCArray Int (Complex Double) {- ^ zx -} ->
   Int {- ^ incx -} ->
   IOCArray Int (Complex Double) {- ^ zy -} ->
   Int {- ^ incy -} ->
   IO ()
swap n zx incx zy incy = do
   zxDim0 <- Call.sizes1 <$> getBounds zx
   zyDim0 <- Call.sizes1 <$> getBounds zy
   Call.assert "swap: 1+(n-1)*abs(incx) == zxDim0" (1+(n-1)*abs(incx) == zxDim0)
   Call.assert "swap: 1+(n-1)*abs(incy) == zyDim0" (1+(n-1)*abs(incy) == zyDim0)
   evalContT $ do
      nPtr <- Call.cint n
      zxPtr <- Call.ioarray zx
      incxPtr <- Call.cint incx
      zyPtr <- Call.ioarray zy
      incyPtr <- Call.cint incy
      liftIO $ FFI.swap nPtr zxPtr incxPtr zyPtr incyPtr

symm ::
   Char {- ^ side -} ->
   Char {- ^ uplo -} ->
   Int {- ^ m -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   CArray (Int,Int) (Complex Double) {- ^ b -} ->
   Complex Double {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ c -} ->
   IO ()
symm side uplo m alpha a b beta c = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   let (bDim0,bDim1) = Call.sizes2 $ bounds b
   (cDim0,cDim1) <- Call.sizes2 <$> getBounds c
   let _ka = aDim0
   let lda = aDim1
   let n = bDim0
   let ldb = bDim1
   let ldc = cDim1
   Call.assert "symm: n == cDim0" (n == cDim0)
   evalContT $ do
      sidePtr <- Call.char side
      uploPtr <- Call.char uplo
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.complexDouble beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.symm sidePtr uploPtr mPtr nPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

syr2k ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Int {- ^ k -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   CArray (Int,Int) (Complex Double) {- ^ b -} ->
   Complex Double {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ c -} ->
   IO ()
syr2k uplo trans k alpha a b beta c = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   let (bDim0,bDim1) = Call.sizes2 $ bounds b
   (cDim0,cDim1) <- Call.sizes2 <$> getBounds c
   let _ka = aDim0
   let lda = aDim1
   let _kb = bDim0
   let ldb = bDim1
   let n = cDim0
   let ldc = cDim1
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.complexDouble beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.syr2k uploPtr transPtr nPtr kPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

syrk ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Int {- ^ k -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   Complex Double {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ c -} ->
   IO ()
syrk uplo trans k alpha a beta c = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   (cDim0,cDim1) <- Call.sizes2 <$> getBounds c
   let _ka = aDim0
   let lda = aDim1
   let n = cDim0
   let ldc = cDim1
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      betaPtr <- Call.complexDouble beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.syrk uploPtr transPtr nPtr kPtr alphaPtr aPtr ldaPtr betaPtr cPtr ldcPtr

tbmv ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Char {- ^ diag -} ->
   Int {- ^ k -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   IOCArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   IO ()
tbmv uplo trans diag k a x incx = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   xDim0 <- Call.sizes1 <$> getBounds x
   let n = aDim0
   let lda = aDim1
   let _xSize = xDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      diagPtr <- Call.char diag
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.ioarray x
      incxPtr <- Call.cint incx
      liftIO $ FFI.tbmv uploPtr transPtr diagPtr nPtr kPtr aPtr ldaPtr xPtr incxPtr

tbsv ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Char {- ^ diag -} ->
   Int {- ^ k -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   IOCArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   IO ()
tbsv uplo trans diag k a x incx = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   xDim0 <- Call.sizes1 <$> getBounds x
   let n = aDim0
   let lda = aDim1
   let _xSize = xDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      diagPtr <- Call.char diag
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.ioarray x
      incxPtr <- Call.cint incx
      liftIO $ FFI.tbsv uploPtr transPtr diagPtr nPtr kPtr aPtr ldaPtr xPtr incxPtr

tpmv ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Char {- ^ diag -} ->
   Int {- ^ n -} ->
   CArray Int (Complex Double) {- ^ ap -} ->
   IOCArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   IO ()
tpmv uplo trans diag n ap x incx = do
   let apDim0 = Call.sizes1 $ bounds ap
   xDim0 <- Call.sizes1 <$> getBounds x
   let _apSize = apDim0
   let _xSize = xDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      diagPtr <- Call.char diag
      nPtr <- Call.cint n
      apPtr <- Call.array ap
      xPtr <- Call.ioarray x
      incxPtr <- Call.cint incx
      liftIO $ FFI.tpmv uploPtr transPtr diagPtr nPtr apPtr xPtr incxPtr

tpsv ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Char {- ^ diag -} ->
   Int {- ^ n -} ->
   CArray Int (Complex Double) {- ^ ap -} ->
   IOCArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   IO ()
tpsv uplo trans diag n ap x incx = do
   let apDim0 = Call.sizes1 $ bounds ap
   xDim0 <- Call.sizes1 <$> getBounds x
   let _apSize = apDim0
   let _xSize = xDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      diagPtr <- Call.char diag
      nPtr <- Call.cint n
      apPtr <- Call.array ap
      xPtr <- Call.ioarray x
      incxPtr <- Call.cint incx
      liftIO $ FFI.tpsv uploPtr transPtr diagPtr nPtr apPtr xPtr incxPtr

trmm ::
   Char {- ^ side -} ->
   Char {- ^ uplo -} ->
   Char {- ^ transa -} ->
   Char {- ^ diag -} ->
   Int {- ^ m -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ b -} ->
   IO ()
trmm side uplo transa diag m alpha a b = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   (bDim0,bDim1) <- Call.sizes2 <$> getBounds b
   let _k = aDim0
   let lda = aDim1
   let n = bDim0
   let ldb = bDim1
   evalContT $ do
      sidePtr <- Call.char side
      uploPtr <- Call.char uplo
      transaPtr <- Call.char transa
      diagPtr <- Call.char diag
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.ioarray b
      ldbPtr <- Call.cint ldb
      liftIO $ FFI.trmm sidePtr uploPtr transaPtr diagPtr mPtr nPtr alphaPtr aPtr ldaPtr bPtr ldbPtr

trmv ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Char {- ^ diag -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   IOCArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   IO ()
trmv uplo trans diag a x incx = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   xDim0 <- Call.sizes1 <$> getBounds x
   let n = aDim0
   let lda = aDim1
   let _xSize = xDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      diagPtr <- Call.char diag
      nPtr <- Call.cint n
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.ioarray x
      incxPtr <- Call.cint incx
      liftIO $ FFI.trmv uploPtr transPtr diagPtr nPtr aPtr ldaPtr xPtr incxPtr

trsm ::
   Char {- ^ side -} ->
   Char {- ^ uplo -} ->
   Char {- ^ transa -} ->
   Char {- ^ diag -} ->
   Int {- ^ m -} ->
   Complex Double {- ^ alpha -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   IOCArray (Int,Int) (Complex Double) {- ^ b -} ->
   IO ()
trsm side uplo transa diag m alpha a b = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   (bDim0,bDim1) <- Call.sizes2 <$> getBounds b
   let _k = aDim0
   let lda = aDim1
   let n = bDim0
   let ldb = bDim1
   evalContT $ do
      sidePtr <- Call.char side
      uploPtr <- Call.char uplo
      transaPtr <- Call.char transa
      diagPtr <- Call.char diag
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.complexDouble alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.ioarray b
      ldbPtr <- Call.cint ldb
      liftIO $ FFI.trsm sidePtr uploPtr transaPtr diagPtr mPtr nPtr alphaPtr aPtr ldaPtr bPtr ldbPtr

trsv ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Char {- ^ diag -} ->
   CArray (Int,Int) (Complex Double) {- ^ a -} ->
   IOCArray Int (Complex Double) {- ^ x -} ->
   Int {- ^ incx -} ->
   IO ()
trsv uplo trans diag a x incx = do
   let (aDim0,aDim1) = Call.sizes2 $ bounds a
   xDim0 <- Call.sizes1 <$> getBounds x
   let n = aDim0
   let lda = aDim1
   let _xSize = xDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      transPtr <- Call.char trans
      diagPtr <- Call.char diag
      nPtr <- Call.cint n
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.ioarray x
      incxPtr <- Call.cint incx
      liftIO $ FFI.trsv uploPtr transPtr diagPtr nPtr aPtr ldaPtr xPtr incxPtr
