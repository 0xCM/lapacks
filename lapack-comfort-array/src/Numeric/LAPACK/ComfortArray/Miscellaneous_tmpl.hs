module Numeric.LAPACK.ComfortArray.Miscellaneous where

import qualified Numeric.LAPACK.FFI.Miscellaneous as FFI
import qualified Numeric.Netlib.ComfortArray.Utility as Call

import Foreign.C.Types (CInt, CChar)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)

