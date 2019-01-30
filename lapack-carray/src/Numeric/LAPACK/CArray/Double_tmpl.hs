module Numeric.LAPACK.CArray.Double where

import qualified Numeric.LAPACK.FFI.Double as FFI
import qualified Numeric.Netlib.CArray.Utility as Call
import Numeric.Netlib.CArray.Utility ((^!))

import Data.Array.IOCArray (IOCArray, getBounds)
import Data.Array.CArray (CArray, bounds)

import Foreign.Storable (peek)
import Foreign.Ptr (Ptr, FunPtr)
import Foreign.C.String (castCCharToChar)
import Foreign.C.Types (CInt)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Applicative (pure, (<*>), (<$>))

