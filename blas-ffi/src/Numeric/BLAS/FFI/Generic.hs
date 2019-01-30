-- Do not edit! Automatically generated by create-lapack-ffi.
module Numeric.BLAS.FFI.Generic (
   axpy,
   copy,
   gbmv,
   gemm,
   gemv,
   gerc,
   geru,
   hbmv,
   hemm,
   hemv,
   her2,
   hpmv,
   hpr2,
   iamax,
   scal,
   swap,
   symm,
   syr2k,
   syrk,
   tbmv,
   tbsv,
   tpmv,
   tpsv,
   trmm,
   trmv,
   trsm,
   trsv,
   ) where

import qualified Numeric.BLAS.FFI.ComplexFloat as C
import qualified Numeric.BLAS.FFI.ComplexDouble as Z
import qualified Numeric.BLAS.FFI.Float as S
import qualified Numeric.BLAS.FFI.Double as D

import qualified Numeric.Netlib.Class as Class

import Foreign.Ptr (Ptr)
import Foreign.C.Types



newtype AXPY a = AXPY {getAXPY :: Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

axpy :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
axpy = getAXPY $ Class.switchFloating (AXPY S.axpy) (AXPY D.axpy) (AXPY C.axpy) (AXPY Z.axpy)


newtype COPY a = COPY {getCOPY :: Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

copy :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
copy = getCOPY $ Class.switchFloating (COPY S.copy) (COPY D.copy) (COPY C.copy) (COPY Z.copy)


newtype GBMV a = GBMV {getGBMV :: Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

gbmv :: Class.Floating a => Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
gbmv = getGBMV $ Class.switchFloating (GBMV S.gbmv) (GBMV D.gbmv) (GBMV C.gbmv) (GBMV Z.gbmv)


newtype GEMM a = GEMM {getGEMM :: Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

gemm :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
gemm = getGEMM $ Class.switchFloating (GEMM S.gemm) (GEMM D.gemm) (GEMM C.gemm) (GEMM Z.gemm)


newtype GEMV a = GEMV {getGEMV :: Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

gemv :: Class.Floating a => Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
gemv = getGEMV $ Class.switchFloating (GEMV S.gemv) (GEMV D.gemv) (GEMV C.gemv) (GEMV Z.gemv)


newtype GERC a = GERC {getGERC :: Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

gerc :: Class.Floating a => Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
gerc = getGERC $ Class.switchFloating (GERC S.ger) (GERC D.ger) (GERC C.gerc) (GERC Z.gerc)


newtype GERU a = GERU {getGERU :: Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

geru :: Class.Floating a => Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
geru = getGERU $ Class.switchFloating (GERU S.ger) (GERU D.ger) (GERU C.geru) (GERU Z.geru)


newtype HBMV a = HBMV {getHBMV :: Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

hbmv :: Class.Floating a => Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
hbmv = getHBMV $ Class.switchFloating (HBMV S.sbmv) (HBMV D.sbmv) (HBMV C.hbmv) (HBMV Z.hbmv)


newtype HEMM a = HEMM {getHEMM :: Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

hemm :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
hemm = getHEMM $ Class.switchFloating (HEMM S.symm) (HEMM D.symm) (HEMM C.hemm) (HEMM Z.hemm)


newtype HEMV a = HEMV {getHEMV :: Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

hemv :: Class.Floating a => Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
hemv = getHEMV $ Class.switchFloating (HEMV S.symv) (HEMV D.symv) (HEMV C.hemv) (HEMV Z.hemv)


newtype HER2 a = HER2 {getHER2 :: Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

her2 :: Class.Floating a => Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
her2 = getHER2 $ Class.switchFloating (HER2 S.syr2) (HER2 D.syr2) (HER2 C.her2) (HER2 Z.her2)


newtype HPMV a = HPMV {getHPMV :: Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

hpmv :: Class.Floating a => Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
hpmv = getHPMV $ Class.switchFloating (HPMV S.spmv) (HPMV D.spmv) (HPMV C.hpmv) (HPMV Z.hpmv)


newtype HPR2 a = HPR2 {getHPR2 :: Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> IO ()}

hpr2 :: Class.Floating a => Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> IO ()
hpr2 = getHPR2 $ Class.switchFloating (HPR2 S.spr2) (HPR2 D.spr2) (HPR2 C.hpr2) (HPR2 Z.hpr2)


newtype IAMAX a = IAMAX {getIAMAX :: Ptr CInt -> Ptr a -> Ptr CInt -> IO CInt}

iamax :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr CInt -> IO CInt
iamax = getIAMAX $ Class.switchFloating (IAMAX S.iamax) (IAMAX D.iamax) (IAMAX C.iamax) (IAMAX Z.iamax)


newtype SCAL a = SCAL {getSCAL :: Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

scal :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
scal = getSCAL $ Class.switchFloating (SCAL S.scal) (SCAL D.scal) (SCAL C.scal) (SCAL Z.scal)


newtype SWAP a = SWAP {getSWAP :: Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

swap :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
swap = getSWAP $ Class.switchFloating (SWAP S.swap) (SWAP D.swap) (SWAP C.swap) (SWAP Z.swap)


newtype SYMM a = SYMM {getSYMM :: Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

symm :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
symm = getSYMM $ Class.switchFloating (SYMM S.symm) (SYMM D.symm) (SYMM C.symm) (SYMM Z.symm)


newtype SYR2K a = SYR2K {getSYR2K :: Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

syr2k :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
syr2k = getSYR2K $ Class.switchFloating (SYR2K S.syr2k) (SYR2K D.syr2k) (SYR2K C.syr2k) (SYR2K Z.syr2k)


newtype SYRK a = SYRK {getSYRK :: Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

syrk :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
syrk = getSYRK $ Class.switchFloating (SYRK S.syrk) (SYRK D.syrk) (SYRK C.syrk) (SYRK Z.syrk)


newtype TBMV a = TBMV {getTBMV :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

tbmv :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
tbmv = getTBMV $ Class.switchFloating (TBMV S.tbmv) (TBMV D.tbmv) (TBMV C.tbmv) (TBMV Z.tbmv)


newtype TBSV a = TBSV {getTBSV :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

tbsv :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
tbsv = getTBSV $ Class.switchFloating (TBSV S.tbsv) (TBSV D.tbsv) (TBSV C.tbsv) (TBSV Z.tbsv)


newtype TPMV a = TPMV {getTPMV :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

tpmv :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
tpmv = getTPMV $ Class.switchFloating (TPMV S.tpmv) (TPMV D.tpmv) (TPMV C.tpmv) (TPMV Z.tpmv)


newtype TPSV a = TPSV {getTPSV :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()}

tpsv :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
tpsv = getTPSV $ Class.switchFloating (TPSV S.tpsv) (TPSV D.tpsv) (TPSV C.tpsv) (TPSV Z.tpsv)


newtype TRMM a = TRMM {getTRMM :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

trmm :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
trmm = getTRMM $ Class.switchFloating (TRMM S.trmm) (TRMM D.trmm) (TRMM C.trmm) (TRMM Z.trmm)


newtype TRMV a = TRMV {getTRMV :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

trmv :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
trmv = getTRMV $ Class.switchFloating (TRMV S.trmv) (TRMV D.trmv) (TRMV C.trmv) (TRMV Z.trmv)


newtype TRSM a = TRSM {getTRSM :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

trsm :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
trsm = getTRSM $ Class.switchFloating (TRSM S.trsm) (TRSM D.trsm) (TRSM C.trsm) (TRSM Z.trsm)


newtype TRSV a = TRSV {getTRSV :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()}

trsv :: Class.Floating a => Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> IO ()
trsv = getTRSV $ Class.switchFloating (TRSV S.trsv) (TRSV D.trsv) (TRSV C.trsv) (TRSV Z.trsv)
