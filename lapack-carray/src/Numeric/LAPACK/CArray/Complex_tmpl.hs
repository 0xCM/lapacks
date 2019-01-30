module Numeric.LAPACK.CArray.Complex where

import qualified Numeric.LAPACK.CArray.ComplexDouble as Z
import qualified Numeric.LAPACK.CArray.ComplexFloat as C
import qualified Numeric.Netlib.Class as Class

import Data.Complex (Complex)

import Data.Array.IOCArray (IOCArray)
import Data.Array.CArray (CArray)

import Foreign.Ptr (Ptr, FunPtr)
import Foreign.C.Types (CInt)

