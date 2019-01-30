module Numeric.LAPACK.CArray.Real where

import qualified Numeric.LAPACK.CArray.Double as D
import qualified Numeric.LAPACK.CArray.Float as S
import qualified Numeric.Netlib.Class as Class

import Data.Array.IOCArray (IOCArray)
import Data.Array.CArray (CArray)

import Foreign.Ptr (Ptr, FunPtr)
import Foreign.C.Types (CInt)

