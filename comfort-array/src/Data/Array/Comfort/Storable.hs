module Data.Array.Comfort.Storable (
   Array,
   shape,
   reshape,
   mapShape,

   (!),
   Array.toList,
   Array.fromList,
   Array.vectorFromList,

   Array.map,
   Array.mapWithIndex,
   (//),
   accumulate,
   fromAssociations,
   ) where

import qualified Data.Array.Comfort.Storable.Mutable.Unchecked as MutArrayNC
import qualified Data.Array.Comfort.Storable.Mutable as MutArray
import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array)

import Foreign.Storable (Storable)

import Control.Monad.ST (runST)

import Data.Foldable (forM_)

import Text.Printf (printf)

import Prelude hiding (map)


shape :: Array sh a -> sh
shape = Array.shape

reshape :: (Shape.C sh0, Shape.C sh1) => sh1 -> Array sh0 a -> Array sh1 a
reshape sh1 arr =
   let n0 = Shape.size $ shape arr
       n1 = Shape.size sh1
   in if n0 == n1
         then Array.reshape sh1 arr
         else error $
              printf
                 ("Array.Comfort.Storable.reshape: " ++
                  "different sizes of old (%d) and new (%d) shape")
                 n0 n1

mapShape ::
   (Shape.C sh0, Shape.C sh1) => (sh0 -> sh1) -> Array sh0 a -> Array sh1 a
mapShape f arr = reshape (f $ shape arr) arr


infixl 9 !

(!) :: (Shape.Indexed sh, Storable a) => Array sh a -> Shape.Index sh -> a
(!) arr ix = runST (do
   marr <- MutArrayNC.unsafeThaw arr
   MutArray.read marr ix)


(//) ::
   (Shape.Indexed sh, Storable a) =>
   Array sh a -> [(Shape.Index sh, a)] -> Array sh a
(//) arr xs = runST (do
   marr <- MutArray.thaw arr
   forM_ xs $ uncurry $ MutArray.write marr
   MutArrayNC.unsafeFreeze marr)

accumulate ::
   (Shape.Indexed sh, Storable a) =>
   (a -> b -> a) -> Array sh a -> [(Shape.Index sh, b)] -> Array sh a
accumulate f arr xs = runST (do
   marr <- MutArray.thaw arr
   forM_ xs $ \(ix,b) -> MutArray.update marr ix $ flip f b
   MutArrayNC.unsafeFreeze marr)

fromAssociations ::
   (Shape.Indexed sh, Storable a) =>
   sh -> a -> [(Shape.Index sh, a)] -> Array sh a
fromAssociations sh a xs = runST (do
   marr <- MutArray.new sh a
   forM_ xs $ uncurry $ MutArray.write marr
   MutArrayNC.unsafeFreeze marr)
