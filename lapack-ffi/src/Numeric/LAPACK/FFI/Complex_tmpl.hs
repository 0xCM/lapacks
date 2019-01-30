module Numeric.LAPACK.FFI.Complex where

import qualified Numeric.LAPACK.FFI.ComplexFloat as C
import qualified Numeric.LAPACK.FFI.ComplexDouble as Z

import qualified Numeric.Netlib.Class as Class

import Data.Complex (Complex)
import Foreign.Ptr (FunPtr, Ptr)
import Foreign.C.Types

