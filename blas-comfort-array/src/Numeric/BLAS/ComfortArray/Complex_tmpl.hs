module Numeric.BLAS.ComfortArray.Complex where

import qualified Numeric.BLAS.ComfortArray.ComplexDouble as Z
import qualified Numeric.BLAS.ComfortArray.ComplexFloat as C
import qualified Numeric.Netlib.Class as Class
import Numeric.Netlib.ComfortArray.Utility (ZeroInt)

import Data.Complex (Complex)

import Data.Array.Comfort.Storable.Mutable (IOArray)
import Data.Array.Comfort.Storable (Array)

import Foreign.C.Types (CInt)

