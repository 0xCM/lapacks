module Numeric.BLAS.ComfortArray.Generic where

import qualified Numeric.BLAS.ComfortArray.ComplexDouble as Z
import qualified Numeric.BLAS.ComfortArray.ComplexFloat as C
import qualified Numeric.BLAS.ComfortArray.Double as D
import qualified Numeric.BLAS.ComfortArray.Float as S
import qualified Numeric.Netlib.Class as Class
import Numeric.Netlib.ComfortArray.Utility (ZeroInt)

import Data.Array.Comfort.Storable.Mutable (IOArray)
import Data.Array.Comfort.Storable (Array)

import Foreign.C.Types (CInt)

