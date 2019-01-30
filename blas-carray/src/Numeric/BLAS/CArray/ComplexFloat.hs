-- Do not edit! Automatically generated by create-lapack-ffi.
module Numeric.BLAS.CArray.ComplexFloat where

import qualified Numeric.BLAS.FFI.ComplexFloat as FFI
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
   Complex Float {- ^ ca -} ->
   CArray Int (Complex Float) {- ^ cx -} ->
   Int {- ^ incx -} ->
   IOCArray Int (Complex Float) {- ^ cy -} ->
   Int {- ^ incy -} ->
   IO ()
axpy n ca cx incx cy incy = do
   let cxDim0 = Call.sizes1 $ bounds cx
   cyDim0 <- Call.sizes1 <$> getBounds cy
   Call.assert "axpy: 1+(n-1)*abs(incx) == cxDim0" (1+(n-1)*abs(incx) == cxDim0)
   Call.assert "axpy: 1+(n-1)*abs(incy) == cyDim0" (1+(n-1)*abs(incy) == cyDim0)
   evalContT $ do
      nPtr <- Call.cint n
      caPtr <- Call.complexFloat ca
      cxPtr <- Call.array cx
      incxPtr <- Call.cint incx
      cyPtr <- Call.ioarray cy
      incyPtr <- Call.cint incy
      liftIO $ FFI.axpy nPtr caPtr cxPtr incxPtr cyPtr incyPtr

casum ::
   Int {- ^ n -} ->
   IOCArray Int (Complex Float) {- ^ cx -} ->
   Int {- ^ incx -} ->
   IO Float
casum n cx incx = do
   cxDim0 <- Call.sizes1 <$> getBounds cx
   Call.assert "casum: 1+(n-1)*abs(incx) == cxDim0" (1+(n-1)*abs(incx) == cxDim0)
   evalContT $ do
      nPtr <- Call.cint n
      cxPtr <- Call.ioarray cx
      incxPtr <- Call.cint incx
      liftIO $ FFI.casum nPtr cxPtr incxPtr

cnrm2 ::
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   IO Float
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
   CArray Int (Complex Float) {- ^ cx -} ->
   Int {- ^ incx -} ->
   Int {- ^ incy -} ->
   IO (CArray Int (Complex Float))
copy n cx incx incy = do
   let cxDim0 = Call.sizes1 $ bounds cx
   Call.assert "copy: 1+(n-1)*abs(incx) == cxDim0" (1+(n-1)*abs(incx) == cxDim0)
   cy <- Call.newArray1 (1+(n-1)*abs(incy))
   evalContT $ do
      nPtr <- Call.cint n
      cxPtr <- Call.array cx
      incxPtr <- Call.cint incx
      cyPtr <- Call.ioarray cy
      incyPtr <- Call.cint incy
      liftIO $ FFI.copy nPtr cxPtr incxPtr cyPtr incyPtr
      liftIO $ Call.freezeArray cy

gbmv ::
   Char {- ^ trans -} ->
   Int {- ^ m -} ->
   Int {- ^ kl -} ->
   Int {- ^ ku -} ->
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Float {- ^ beta -} ->
   IOCArray Int (Complex Float) {- ^ y -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexFloat beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.gbmv transPtr mPtr nPtr klPtr kuPtr alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

gemm ::
   Char {- ^ transa -} ->
   Char {- ^ transb -} ->
   Int {- ^ m -} ->
   Int {- ^ k -} ->
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   CArray (Int,Int) (Complex Float) {- ^ b -} ->
   Complex Float {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ c -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.complexFloat beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.gemm transaPtr transbPtr mPtr nPtr kPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

gemv ::
   Char {- ^ trans -} ->
   Int {- ^ m -} ->
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Float {- ^ beta -} ->
   IOCArray Int (Complex Float) {- ^ y -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexFloat beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.gemv transPtr mPtr nPtr alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

gerc ::
   Int {- ^ m -} ->
   Complex Float {- ^ alpha -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   CArray Int (Complex Float) {- ^ y -} ->
   Int {- ^ incy -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ a -} ->
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
      alphaPtr <- Call.complexFloat alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      yPtr <- Call.array y
      incyPtr <- Call.cint incy
      aPtr <- Call.ioarray a
      ldaPtr <- Call.cint lda
      liftIO $ FFI.gerc mPtr nPtr alphaPtr xPtr incxPtr yPtr incyPtr aPtr ldaPtr

geru ::
   Int {- ^ m -} ->
   Complex Float {- ^ alpha -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   CArray Int (Complex Float) {- ^ y -} ->
   Int {- ^ incy -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ a -} ->
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
      alphaPtr <- Call.complexFloat alpha
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
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Float {- ^ beta -} ->
   IOCArray Int (Complex Float) {- ^ y -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexFloat beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.hbmv uploPtr nPtr kPtr alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

hemm ::
   Char {- ^ side -} ->
   Char {- ^ uplo -} ->
   Int {- ^ m -} ->
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   CArray (Int,Int) (Complex Float) {- ^ b -} ->
   Complex Float {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ c -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.complexFloat beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.hemm sidePtr uploPtr mPtr nPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

hemv ::
   Char {- ^ uplo -} ->
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Float {- ^ beta -} ->
   IOCArray Int (Complex Float) {- ^ y -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexFloat beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.hemv uploPtr nPtr alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

her ::
   Char {- ^ uplo -} ->
   Float {- ^ alpha -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ a -} ->
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
      alphaPtr <- Call.float alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      aPtr <- Call.ioarray a
      ldaPtr <- Call.cint lda
      liftIO $ FFI.her uploPtr nPtr alphaPtr xPtr incxPtr aPtr ldaPtr

her2 ::
   Char {- ^ uplo -} ->
   Complex Float {- ^ alpha -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   CArray Int (Complex Float) {- ^ y -} ->
   Int {- ^ incy -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ a -} ->
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
      alphaPtr <- Call.complexFloat alpha
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
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   CArray (Int,Int) (Complex Float) {- ^ b -} ->
   Float {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ c -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.float beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.her2k uploPtr transPtr nPtr kPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

herk ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Int {- ^ k -} ->
   Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   Float {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ c -} ->
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
      alphaPtr <- Call.float alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      betaPtr <- Call.float beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.herk uploPtr transPtr nPtr kPtr alphaPtr aPtr ldaPtr betaPtr cPtr ldcPtr

hpmv ::
   Char {- ^ uplo -} ->
   Int {- ^ n -} ->
   Complex Float {- ^ alpha -} ->
   CArray Int (Complex Float) {- ^ ap -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   Complex Float {- ^ beta -} ->
   IOCArray Int (Complex Float) {- ^ y -} ->
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
      alphaPtr <- Call.complexFloat alpha
      apPtr <- Call.array ap
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      betaPtr <- Call.complexFloat beta
      yPtr <- Call.ioarray y
      incyPtr <- Call.cint incy
      liftIO $ FFI.hpmv uploPtr nPtr alphaPtr apPtr xPtr incxPtr betaPtr yPtr incyPtr

hpr ::
   Char {- ^ uplo -} ->
   Int {- ^ n -} ->
   Float {- ^ alpha -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   IOCArray Int (Complex Float) {- ^ ap -} ->
   IO ()
hpr uplo n alpha x incx ap = do
   let xDim0 = Call.sizes1 $ bounds x
   apDim0 <- Call.sizes1 <$> getBounds ap
   let _xSize = xDim0
   let _apSize = apDim0
   evalContT $ do
      uploPtr <- Call.char uplo
      nPtr <- Call.cint n
      alphaPtr <- Call.float alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      apPtr <- Call.ioarray ap
      liftIO $ FFI.hpr uploPtr nPtr alphaPtr xPtr incxPtr apPtr

hpr2 ::
   Char {- ^ uplo -} ->
   Int {- ^ n -} ->
   Complex Float {- ^ alpha -} ->
   CArray Int (Complex Float) {- ^ x -} ->
   Int {- ^ incx -} ->
   CArray Int (Complex Float) {- ^ y -} ->
   Int {- ^ incy -} ->
   IOCArray Int (Complex Float) {- ^ ap -} ->
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
      alphaPtr <- Call.complexFloat alpha
      xPtr <- Call.array x
      incxPtr <- Call.cint incx
      yPtr <- Call.array y
      incyPtr <- Call.cint incy
      apPtr <- Call.ioarray ap
      liftIO $ FFI.hpr2 uploPtr nPtr alphaPtr xPtr incxPtr yPtr incyPtr apPtr

iamax ::
   Int {- ^ n -} ->
   CArray Int (Complex Float) {- ^ cx -} ->
   Int {- ^ incx -} ->
   IO CInt
iamax n cx incx = do
   let cxDim0 = Call.sizes1 $ bounds cx
   Call.assert "iamax: 1+(n-1)*abs(incx) == cxDim0" (1+(n-1)*abs(incx) == cxDim0)
   evalContT $ do
      nPtr <- Call.cint n
      cxPtr <- Call.array cx
      incxPtr <- Call.cint incx
      liftIO $ FFI.iamax nPtr cxPtr incxPtr

rotg ::
   Complex Float {- ^ ca -} ->
   Complex Float {- ^ cb -} ->
   IO (Float, Complex Float)
rotg ca cb = do
   evalContT $ do
      caPtr <- Call.complexFloat ca
      cbPtr <- Call.complexFloat cb
      cPtr <- Call.alloca
      sPtr <- Call.alloca
      liftIO $ FFI.rotg caPtr cbPtr cPtr sPtr
      liftIO $ pure (,)
         <*> peek cPtr
         <*> peek sPtr

rrot ::
   Int {- ^ n -} ->
   IOCArray Int (Complex Float) {- ^ cx -} ->
   Int {- ^ incx -} ->
   IOCArray Int (Complex Float) {- ^ cy -} ->
   Int {- ^ incy -} ->
   Float {- ^ c -} ->
   Float {- ^ s -} ->
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
      cPtr <- Call.float c
      sPtr <- Call.float s
      liftIO $ FFI.rrot nPtr cxPtr incxPtr cyPtr incyPtr cPtr sPtr

rscal ::
   Int {- ^ n -} ->
   Float {- ^ sa -} ->
   IOCArray Int (Complex Float) {- ^ cx -} ->
   Int {- ^ incx -} ->
   IO ()
rscal n sa cx incx = do
   cxDim0 <- Call.sizes1 <$> getBounds cx
   Call.assert "rscal: 1+(n-1)*abs(incx) == cxDim0" (1+(n-1)*abs(incx) == cxDim0)
   evalContT $ do
      nPtr <- Call.cint n
      saPtr <- Call.float sa
      cxPtr <- Call.ioarray cx
      incxPtr <- Call.cint incx
      liftIO $ FFI.rscal nPtr saPtr cxPtr incxPtr

scal ::
   Int {- ^ n -} ->
   Complex Float {- ^ ca -} ->
   IOCArray Int (Complex Float) {- ^ cx -} ->
   Int {- ^ incx -} ->
   IO ()
scal n ca cx incx = do
   cxDim0 <- Call.sizes1 <$> getBounds cx
   Call.assert "scal: 1+(n-1)*abs(incx) == cxDim0" (1+(n-1)*abs(incx) == cxDim0)
   evalContT $ do
      nPtr <- Call.cint n
      caPtr <- Call.complexFloat ca
      cxPtr <- Call.ioarray cx
      incxPtr <- Call.cint incx
      liftIO $ FFI.scal nPtr caPtr cxPtr incxPtr

swap ::
   Int {- ^ n -} ->
   IOCArray Int (Complex Float) {- ^ cx -} ->
   Int {- ^ incx -} ->
   IOCArray Int (Complex Float) {- ^ cy -} ->
   Int {- ^ incy -} ->
   IO ()
swap n cx incx cy incy = do
   cxDim0 <- Call.sizes1 <$> getBounds cx
   cyDim0 <- Call.sizes1 <$> getBounds cy
   Call.assert "swap: 1+(n-1)*abs(incx) == cxDim0" (1+(n-1)*abs(incx) == cxDim0)
   Call.assert "swap: 1+(n-1)*abs(incy) == cyDim0" (1+(n-1)*abs(incy) == cyDim0)
   evalContT $ do
      nPtr <- Call.cint n
      cxPtr <- Call.ioarray cx
      incxPtr <- Call.cint incx
      cyPtr <- Call.ioarray cy
      incyPtr <- Call.cint incy
      liftIO $ FFI.swap nPtr cxPtr incxPtr cyPtr incyPtr

symm ::
   Char {- ^ side -} ->
   Char {- ^ uplo -} ->
   Int {- ^ m -} ->
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   CArray (Int,Int) (Complex Float) {- ^ b -} ->
   Complex Float {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ c -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.complexFloat beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.symm sidePtr uploPtr mPtr nPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

syr2k ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Int {- ^ k -} ->
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   CArray (Int,Int) (Complex Float) {- ^ b -} ->
   Complex Float {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ c -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.array b
      ldbPtr <- Call.cint ldb
      betaPtr <- Call.complexFloat beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.syr2k uploPtr transPtr nPtr kPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr

syrk ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Int {- ^ k -} ->
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   Complex Float {- ^ beta -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ c -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      betaPtr <- Call.complexFloat beta
      cPtr <- Call.ioarray c
      ldcPtr <- Call.cint ldc
      liftIO $ FFI.syrk uploPtr transPtr nPtr kPtr alphaPtr aPtr ldaPtr betaPtr cPtr ldcPtr

tbmv ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Char {- ^ diag -} ->
   Int {- ^ k -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   IOCArray Int (Complex Float) {- ^ x -} ->
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
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   IOCArray Int (Complex Float) {- ^ x -} ->
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
   CArray Int (Complex Float) {- ^ ap -} ->
   IOCArray Int (Complex Float) {- ^ x -} ->
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
   CArray Int (Complex Float) {- ^ ap -} ->
   IOCArray Int (Complex Float) {- ^ x -} ->
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
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ b -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.ioarray b
      ldbPtr <- Call.cint ldb
      liftIO $ FFI.trmm sidePtr uploPtr transaPtr diagPtr mPtr nPtr alphaPtr aPtr ldaPtr bPtr ldbPtr

trmv ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Char {- ^ diag -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   IOCArray Int (Complex Float) {- ^ x -} ->
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
   Complex Float {- ^ alpha -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   IOCArray (Int,Int) (Complex Float) {- ^ b -} ->
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
      alphaPtr <- Call.complexFloat alpha
      aPtr <- Call.array a
      ldaPtr <- Call.cint lda
      bPtr <- Call.ioarray b
      ldbPtr <- Call.cint ldb
      liftIO $ FFI.trsm sidePtr uploPtr transaPtr diagPtr mPtr nPtr alphaPtr aPtr ldaPtr bPtr ldbPtr

trsv ::
   Char {- ^ uplo -} ->
   Char {- ^ trans -} ->
   Char {- ^ diag -} ->
   CArray (Int,Int) (Complex Float) {- ^ a -} ->
   IOCArray Int (Complex Float) {- ^ x -} ->
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