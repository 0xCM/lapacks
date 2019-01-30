module Numeric.BLAS.CArray.ComplexFloat where

import qualified Numeric.BLAS.FFI.ComplexFloat as FFI
import qualified Numeric.Netlib.CArray.Utility as Call

import Data.Array.IOCArray (IOCArray, getBounds)
import Data.Array.CArray (CArray, bounds)

import Data.Complex (Complex)

import Foreign.Storable.Complex ()
import Foreign.Storable (peek)
import Foreign.C.Types (CInt)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Applicative (pure, (<$>), (<*>))

