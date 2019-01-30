{-# LANGUAGE TypeFamilies #-}
module Data.RefTF where

import Data.IORef (IORef, newIORef, readIORef, writeIORef, )
import Data.STRef (STRef, newSTRef, readSTRef, writeSTRef, )
import Control.Concurrent.STM.TVar (TVar, newTVar, readTVar, writeTVar, )
import Control.Concurrent.STM (STM, )
import Control.Monad.ST (ST)
import Prelude hiding (read)


modify :: C m => Ref m a -> (a -> a) -> m ()
modify ref f = write ref . f =<< read ref


class Monad m => C m where
   type Ref m :: * -> *
   new :: a -> m (Ref m a)
   write :: Ref m a -> a -> m ()
   read :: Ref m a -> m a

instance C IO where
   type Ref IO = IORef
   new = newIORef
   write = writeIORef
   read = readIORef

instance C (ST s) where
   type Ref (ST s) = STRef s
   new = newSTRef
   write = writeSTRef
   read = readSTRef

instance C STM where
   type Ref STM = TVar
   new = newTVar
   write = writeTVar
   read = readTVar
