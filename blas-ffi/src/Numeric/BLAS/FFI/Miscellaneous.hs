-- Do not edit! Automatically generated by create-lapack-ffi.
{-# LANGUAGE ForeignFunctionInterface #-}
module Numeric.BLAS.FFI.Miscellaneous where

import Foreign.Ptr (Ptr)
import Foreign.C.Types


foreign import ccall "lsame_"
   lsame :: Ptr CChar -> Ptr CChar -> IO Bool
