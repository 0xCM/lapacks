module Numeric.LAPACK.ComfortArray.ComplexDouble where

import qualified Numeric.LAPACK.FFI.ComplexDouble as FFI
import qualified Numeric.Netlib.ComfortArray.Utility as Call
import Numeric.Netlib.ComfortArray.Utility (ZeroInt, (^!))

import qualified Data.Array.Comfort.Storable.Mutable as MutArray
import qualified Data.Array.Comfort.Storable as Array
import Data.Array.Comfort.Storable.Mutable (IOArray)
import Data.Array.Comfort.Storable (Array)

import Data.Complex (Complex)

import Foreign.Storable.Complex ()
import Foreign.Storable (peek)
import Foreign.Ptr (Ptr, FunPtr)
import Foreign.C.String (castCCharToChar)
import Foreign.C.Types (CInt)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Applicative (pure, (<*>))

