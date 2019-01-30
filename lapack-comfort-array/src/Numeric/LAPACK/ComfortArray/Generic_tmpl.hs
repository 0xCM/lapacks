module Numeric.LAPACK.ComfortArray.Generic where

import qualified Numeric.LAPACK.ComfortArray.ComplexDouble as Z
import qualified Numeric.LAPACK.ComfortArray.ComplexFloat as C
import qualified Numeric.LAPACK.ComfortArray.Double as D
import qualified Numeric.LAPACK.ComfortArray.Float as S
import qualified Numeric.Netlib.Class as Class
import Numeric.Netlib.ComfortArray.Utility (ZeroInt)

import Data.Array.Comfort.Storable.Mutable (IOArray)
import Data.Array.Comfort.Storable (Array)

import Foreign.C.Types (CInt)

