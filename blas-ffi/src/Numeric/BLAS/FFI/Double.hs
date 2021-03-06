-- Do not edit! Automatically generated by create-lapack-ffi.
{-# LANGUAGE ForeignFunctionInterface #-}
module Numeric.BLAS.FFI.Double where

import Foreign.Ptr (Ptr)
import Foreign.C.Types


foreign import ccall "dasum_"
   asum :: Ptr CInt -> Ptr Double -> Ptr CInt -> IO Double

foreign import ccall "daxpy_"
   axpy :: Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dcopy_"
   copy :: Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "ddot_"
   dot :: Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO Double

foreign import ccall "dgbmv_"
   gbmv :: Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dgemm_"
   gemm :: Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dgemv_"
   gemv :: Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dger_"
   ger :: Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dsbmv_"
   sbmv :: Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dsymv_"
   symv :: Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dsyr_"
   syr :: Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dsyr2_"
   syr2 :: Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dspmv_"
   spmv :: Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dspr_"
   spr :: Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> IO ()

foreign import ccall "dspr2_"
   spr2 :: Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> IO ()

foreign import ccall "idamax_"
   iamax :: Ptr CInt -> Ptr Double -> Ptr CInt -> IO CInt

foreign import ccall "dnrm2_"
   nrm2 :: Ptr CInt -> Ptr Double -> Ptr CInt -> IO Double

foreign import ccall "drot_"
   rot :: Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall "drotg_"
   rotg :: Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall "drotm_"
   rotm :: Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> IO ()

foreign import ccall "drotmg_"
   rotmg :: Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall "dscal_"
   scal :: Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dsdot_"
   sdot :: Ptr CInt -> Ptr Float -> Ptr CInt -> Ptr Float -> Ptr CInt -> IO Double

foreign import ccall "dswap_"
   swap :: Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dsymm_"
   symm :: Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dsyr2k_"
   syr2k :: Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dsyrk_"
   syrk :: Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dtbmv_"
   tbmv :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dtbsv_"
   tbsv :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dtpmv_"
   tpmv :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dtpsv_"
   tpsv :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dtrmm_"
   trmm :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dtrmv_"
   trmv :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dtrsm_"
   trsm :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr Double -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()

foreign import ccall "dtrsv_"
   trsv :: Ptr CChar -> Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr Double -> Ptr CInt -> Ptr Double -> Ptr CInt -> IO ()
