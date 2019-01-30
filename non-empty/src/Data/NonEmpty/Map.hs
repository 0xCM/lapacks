module Data.NonEmpty.Map (
   T,
   insert,
   singleton,
   member,
   size,
   elems,
   keys,
   keysSet,
   lookup,
   minViewWithKey,
   maxViewWithKey,
   fromList,
   toAscList,
   fetch,
   flatten,
   union,
   unionLeft,
   unionRight,
   map,
   mapWithKey,
   ) where

import qualified Data.NonEmpty.Set as NonEmptySet
import qualified Data.NonEmpty.Class as C
import qualified Data.NonEmpty as NonEmpty

import qualified Data.Map as Map
import Data.Map (Map, )

import Control.Monad (mzero, )
import Control.Applicative (liftA2, )
import Control.DeepSeq (NFData, rnf, )
import Data.Traversable (Traversable, traverse, )
import Data.Foldable (Foldable, foldMap, )
import Data.Monoid (mappend, )
import Data.Maybe (fromMaybe, )
import Data.Tuple.HT (forcePair, mapSnd, )
import Data.Ord.HT (comparing, )

import Prelude hiding (map, lookup, )


{-
The first field will always contain the smallest element.
-}
data T k a = Cons (k, a) (Map k a)
   deriving (Eq, Ord)

instance (Show k, Show a) => Show (T k a) where
   showsPrec p xs =
      showParen (p>10) $
         showString "NonEmptyMap.fromList " .
         showsPrec 11 (toAscList xs)

instance (Ord k) => Functor (T k) where
   fmap = map

instance (Ord k) => Foldable (T k) where
   foldMap f (Cons x xs) = mappend (f (snd x)) (foldMap f xs)

instance (Ord k) => Traversable (T k) where
   traverse f (Cons x xs) =
      liftA2 Cons (fmap ((,) (fst x)) $ f (snd x)) (traverse f xs)

instance (NFData k, NFData a) => NFData (T k a) where
   rnf = C.rnf

instance (NFData k) => C.NFData (T k) where
   rnf (Cons x xs) = rnf (x, C.rnf xs)


insert :: Ord k => k -> a -> Map k a -> T k a
insert = curry $ insertGen fst

insertGen :: Ord k => (((k,a),(k,a)) -> (k,a)) -> (k,a) -> Map k a -> T k a
insertGen select y xt =
   uncurry Cons $
   fromMaybe (y, xt) $ do
      (x,xs) <- Map.minViewWithKey xt
      case comparing fst y x of
         GT -> return (x, uncurry Map.insert y xs)
         EQ -> return (select (y,x), xs)
         LT -> mzero

singleton :: k -> a -> T k a
singleton k a = Cons (k,a) Map.empty

member :: (Ord k) => k -> T k a -> Bool
member y (Cons x xs) =
   y == fst x || Map.member y xs

size :: T k a -> Int
size (Cons _ xs) = 1 + Map.size xs

elems :: T k a -> NonEmpty.T [] a
elems (Cons x xs) = NonEmpty.cons (snd x) (Map.elems xs)

keys :: T k a -> NonEmpty.T [] k
keys (Cons x xs) = NonEmpty.cons (fst x) (Map.keys xs)

-- 'insert' could be optimized to 'Cons'
keysSet :: (Ord k) => T k a -> NonEmptySet.T k
keysSet (Cons x xs) = NonEmptySet.insert (fst x) (Map.keysSet xs)

lookup :: (Ord k) => k -> T k a -> Maybe a
lookup y (Cons x xs) =
   if y == fst x
     then Just $ snd x
     else Map.lookup y xs

minViewWithKey :: T k a -> ((k,a), Map k a)
minViewWithKey (Cons x xs) = (x,xs)

maxViewWithKey :: (Ord k) => T k a -> ((k,a), Map k a)
maxViewWithKey (Cons x xs) =
   forcePair $
   case Map.maxViewWithKey xs of
      Nothing -> (x,xs)
      Just (y,ys) -> (y, uncurry Map.insert x ys)

fromList :: (Ord k) => NonEmpty.T [] (k,a) -> T k a
fromList (NonEmpty.Cons x xs) = uncurry insert x $ Map.fromList xs

toAscList :: T k a -> NonEmpty.T [] (k,a)
toAscList (Cons x xs) = NonEmpty.cons x $ Map.toAscList xs

fetch :: (Ord k) => Map k a -> Maybe (T k a)
fetch  =  fmap (uncurry Cons) . Map.minViewWithKey

flatten :: (Ord k) => T k a -> Map k a
flatten (Cons x xs) = uncurry Map.insert x xs

union :: (Ord k) => T k a -> T k a -> T k a
union (Cons x xs) (Cons y ys) =
   uncurry Cons $
   case Map.union xs ys of
      zs ->
         case comparing fst x y of
            LT -> (x, Map.union zs $ uncurry Map.singleton y)
            GT -> (y, uncurry Map.insert x zs)
            EQ -> (x, zs)

unionLeft :: (Ord k) => Map k a -> T k a -> T k a
unionLeft xs (Cons y ys) =
   insertGen snd y $ Map.union xs ys

unionRight :: (Ord k) => T k a -> Map k a -> T k a
unionRight (Cons x xs) ys =
   insertGen fst x $ Map.union xs ys

map :: (Ord k) => (a -> b) -> T k a -> T k b
map f (Cons x xs) = Cons (mapSnd f x) (Map.map f xs)

mapWithKey :: (Ord k) => (k -> a -> b) -> T k a -> T k b
mapWithKey f (Cons x@(k,_a) xs) = Cons (k, uncurry f x) (Map.mapWithKey f xs)
