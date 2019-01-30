module Numeric.BLAS.ComfortArray.Miscellaneous where

import qualified Numeric.BLAS.FFI.Miscellaneous as FFI
import qualified Numeric.Netlib.ComfortArray.Utility as Call

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)

