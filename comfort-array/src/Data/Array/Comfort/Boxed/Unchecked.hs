{-# LANGUAGE TypeFamilies #-}
module Data.Array.Comfort.Boxed.Unchecked where

import qualified Data.Array.Comfort.Shape as Shape
import qualified Data.Primitive.Array as Prim

import qualified Control.Monad.ST.Strict as STStrict
import qualified Control.Monad.Trans.Class as MT
import qualified Control.Monad.Trans.State as MS
import Control.Monad (liftM)
import Control.Applicative ((<$>))
import Control.DeepSeq (NFData, rnf)

import qualified Data.Traversable as Trav
import qualified Data.Foldable as Fold
import qualified Data.List as List
import Prelude hiding (map, )


data Array sh a =
   Array {
      shape :: sh,
      buffer :: Prim.Array a
   }

instance (Shape.C sh, Show sh, Show a) => Show (Array sh a) where
   show arr =
      "BoxedArray.fromList " ++
      showsPrec 11 (shape arr) (' ' : show (toListLazy arr))


instance (Shape.C sh, NFData sh, NFData a) => NFData (Array sh a) where
   rnf a@(Array sh _arr) = rnf (sh, toListLazy a)

instance (Shape.C sh) => Functor (Array sh) where
   fmap = map

instance (Shape.C sh) => Fold.Foldable (Array sh) where
   fold = Fold.fold . buffer
   foldMap f = Fold.foldMap f . buffer
   foldl f a = Fold.foldl f a . buffer
   foldr f a = Fold.foldr f a . buffer
   foldl1 f = Fold.foldl1 f . buffer
   foldr1 f = Fold.foldr1 f . buffer

instance (Shape.C sh) => Trav.Traversable (Array sh) where
   traverse f (Array sh arr) = Array sh <$> Trav.traverse f arr
   sequenceA (Array sh arr) = Array sh <$> Trav.sequenceA arr
   mapM f (Array sh arr) = liftM (Array sh) $ Trav.mapM f arr
   sequence (Array sh arr) = liftM (Array sh) $ Trav.sequence arr


-- add assertion, at least in an exposed version
reshape :: sh1 -> Array sh0 a -> Array sh1 a
reshape sh (Array _ arr) = Array sh arr

mapShape :: (sh0 -> sh1) -> Array sh0 a -> Array sh1 a
mapShape f (Array sh arr) = Array (f sh) arr


infixl 9 !

(!) :: (Shape.Indexed sh) => Array sh a -> Shape.Index sh -> a
(!) (Array sh arr) ix = Prim.indexArray arr $ Shape.uncheckedOffset sh ix

toListLazy :: (Shape.C sh) => Array sh a -> [a]
toListLazy (Array sh arr) =
   List.map (Prim.indexArray arr) $ take (Shape.size sh) [0..]

toList :: (Shape.C sh) => Array sh a -> [a]
toList (Array sh arr) =
   STStrict.runST (mapM (Prim.indexArrayM arr) $ take (Shape.size sh) [0..])

fromList :: (Shape.C sh) => sh -> [a] -> Array sh a
fromList sh xs = Array sh $ Prim.fromListN (Shape.size sh) xs

vectorFromList :: [a] -> Array (Shape.ZeroBased Int) a
vectorFromList xs =
   let arr = Prim.fromList xs
   in Array (Shape.ZeroBased $ Prim.sizeofArray arr) arr

map :: (Shape.C sh) => (a -> b) -> Array sh a -> Array sh b
map f (Array sh arr) = Array sh $ Prim.mapArray' f arr

zipWith ::
   (Shape.C sh) => (a -> b -> c) -> Array sh a -> Array sh b -> Array sh c
zipWith f (Array sha arra) (Array _shb arrb) =
   Array sha $
   STStrict.runST
      (flip MS.evalStateT 0 $
       Prim.traverseArrayP
         (\a -> do
            k <- MS.get
            b <- MT.lift $ Prim.indexArrayM arrb k
            MS.put (k+1)
            return $ f a b)
         arra)
