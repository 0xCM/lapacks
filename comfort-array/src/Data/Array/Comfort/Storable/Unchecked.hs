{-# LANGUAGE TypeFamilies #-}
{- |
The functions in this module miss any bound checking.
-}
module Data.Array.Comfort.Storable.Unchecked (
   Priv.Array(Array, shape, buffer),
   Priv.reshape,
   Priv.mapShape,

   (Priv.!),
   unsafeCreate,
   unsafeCreateWithSize,
   unsafeCreateWithSizeAndResult,
   Priv.toList,
   Priv.fromList,
   Priv.vectorFromList,

   map,
   mapWithIndex,
   (Priv.//),
   Priv.accumulate,
   Priv.fromAssociations,
   ) where

import qualified Data.Array.Comfort.Storable.Unchecked.Monadic as Monadic
import qualified Data.Array.Comfort.Storable.Private as Priv
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Private (Array(Array))

import Foreign.Marshal.Array (advancePtr)
import Foreign.Storable (Storable, poke, peek)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Ptr (Ptr)

import Control.Monad.ST (runST)

import Prelude hiding (map)


unsafeCreate ::
   (Shape.C sh, Storable a) =>
   sh -> (Ptr a -> IO ()) -> Array sh a
unsafeCreate sh arr = runST (Monadic.unsafeCreate sh arr)

unsafeCreateWithSize ::
   (Shape.C sh, Storable a) =>
   sh -> (Int -> Ptr a -> IO ()) -> Array sh a
unsafeCreateWithSize sh arr = runST (Monadic.unsafeCreateWithSize sh arr)

unsafeCreateWithSizeAndResult ::
   (Shape.C sh, Storable a) =>
   sh -> (Int -> Ptr a -> IO b) -> (Array sh a, b)
unsafeCreateWithSizeAndResult sh arr =
   runST (Monadic.unsafeCreateWithSizeAndResult sh arr)


map ::
   (Shape.C sh, Storable a, Storable b) =>
   (a -> b) -> Array sh a -> Array sh b
map f (Array sh a) =
   unsafeCreate sh $ \dstPtr ->
   withForeignPtr a $ \srcPtr ->
   sequence_ $ take (Shape.size sh) $
      zipWith
         (\src dst -> poke dst . f =<< peek src)
         (iterate (flip advancePtr 1) srcPtr)
         (iterate (flip advancePtr 1) dstPtr)

mapWithIndex ::
   (Shape.Indexed sh, Shape.Index sh ~ ix, Storable a, Storable b) =>
   (ix -> a -> b) -> Array sh a -> Array sh b
mapWithIndex f (Array sh a) =
   unsafeCreate sh $ \dstPtr ->
   withForeignPtr a $ \srcPtr ->
   sequence_ $
      zipWith3
         (\ix src dst -> poke dst . f ix =<< peek src)
         (Shape.indices sh)
         (iterate (flip advancePtr 1) srcPtr)
         (iterate (flip advancePtr 1) dstPtr)
