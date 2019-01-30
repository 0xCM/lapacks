module Data.Array.Comfort.Boxed (
   Array,
   shape,
   (!),
   Array.toList,
   toAssociations,
   Array.fromList,
   Array.vectorFromList,

   Array.map,
   zipWith,
   (//),
   accumulate,
   fromAssociations,
   ) where

import qualified Data.Array.Comfort.Boxed.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Boxed.Unchecked (Array(Array))

import qualified Data.Primitive.Array as Prim

import qualified Control.Monad.Primitive as PrimM
import Control.Monad.ST (runST)
import Control.Applicative ((<$>))

import Data.Foldable (forM_)

import Prelude hiding (zipWith)


shape :: Array.Array sh a -> sh
shape = Array.shape


infixl 9 !

(!) :: (Shape.Indexed sh) => Array sh a -> Shape.Index sh -> a
(!) (Array sh arr) ix =
   if Shape.inBounds sh ix
      then Prim.indexArray arr $ Shape.offset sh ix
      else error "Array.Comfort.Boxed.!: index out of bounds"


zipWith ::
   (Shape.C sh, Eq sh) =>
   (a -> b -> c) -> Array sh a -> Array sh b -> Array sh c
zipWith f a b =
   if shape a == shape b
      then Array.zipWith f a b
      else error "zipWith: shapes mismatch"


(//) ::
   (Shape.Indexed sh) => Array sh a -> [(Shape.Index sh, a)] -> Array sh a
(//) (Array sh arr) xs = runST (do
   marr <- Prim.thawArray arr 0 (Shape.size sh)
   forM_ xs $ \(ix,a) -> Prim.writeArray marr (Shape.offset sh ix) a
   Array sh <$> Prim.unsafeFreezeArray marr)

accumulate ::
   (Shape.Indexed sh) =>
   (a -> b -> a) -> Array sh a -> [(Shape.Index sh, b)] -> Array sh a
accumulate f (Array sh arr) xs = runST (do
   marr <- Prim.thawArray arr 0 (Shape.size sh)
   forM_ xs $ \(ix,b) -> updateArray marr (Shape.offset sh ix) $ flip f b
   Array sh <$> Prim.unsafeFreezeArray marr)

updateArray ::
   PrimM.PrimMonad m =>
   Prim.MutableArray (PrimM.PrimState m) a -> Int -> (a -> a) -> m ()
updateArray marr k f = Prim.writeArray marr k . f =<< Prim.readArray marr k

toAssociations :: (Shape.Indexed sh) => Array sh a -> [(Shape.Index sh, a)]
toAssociations arr = zip (Shape.indices $ shape arr) (Array.toList arr)

fromAssociations ::
   (Shape.Indexed sh) => sh -> a -> [(Shape.Index sh, a)] -> Array sh a
fromAssociations sh a xs = runST (do
   marr <- Prim.newArray (Shape.size sh) a
   forM_ xs $ \(ix,x) -> Prim.writeArray marr (Shape.offset sh ix) x
   Array sh <$> Prim.unsafeFreezeArray marr)
