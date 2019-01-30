module Numeric.LAPACK.ComfortArray.Complex where

import qualified Numeric.LAPACK.ComfortArray.ComplexDouble as Z
import qualified Numeric.LAPACK.ComfortArray.ComplexFloat as C
import qualified Numeric.Netlib.Class as Class
import Numeric.Netlib.ComfortArray.Utility (ZeroInt)

import Data.Complex (Complex)

import Data.Array.Comfort.Storable.Mutable (IOArray)
import Data.Array.Comfort.Storable (Array)

import Foreign.Ptr (Ptr, FunPtr)
import Foreign.C.Types (CInt)

