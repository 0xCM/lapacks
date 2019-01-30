{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Data.Array.Comfort.Shape (
   C(..),
   Indexed(..),
   InvIndexed(..),

   ZeroBased(..),
   OneBased(..),

   Range(..),
   Shifted(..),
   Enumeration(..),
   Deferred(..), DeferredIndex, deferIndex, revealIndex,

   (:+:)(..),

   Triangular(..), Lower(Lower), Upper(Upper),
   LowerTriangular, UpperTriangular,
   lowerTriangular, upperTriangular,
   triangleSize, triangleRoot,
   ) where

import qualified Foreign.Storable.Newtype as Store
import Foreign.Storable
         (Storable, sizeOf, alignment, poke, peek, pokeElemOff, peekElemOff)
import Foreign.Ptr (Ptr, castPtr)

import qualified GHC.Arr as Ix

import qualified Control.Monad.Trans.State as MS
import qualified Control.Monad.HT as Monad
import qualified Control.Applicative.Backwards as Back
import Control.DeepSeq (NFData, rnf)
import Control.Applicative (Applicative, pure, liftA2, liftA3, (<*>))
import Control.Applicative (Const(Const, getConst))

import Text.Printf (printf)

import qualified Data.NonEmpty as NonEmpty
import Data.List.HT (tails)
import Data.Tuple.HT (mapSnd, mapPair, swap, fst3, snd3, thd3)


class C sh where
   -- Ix.rangeSize
   size :: sh -> Int
   -- Ix.unsafeRangeSize
   uncheckedSize :: sh -> Int
   uncheckedSize = size

class C sh => Indexed sh where
   {-# MINIMAL indices, (sizeOffset|offset), inBounds #-}
   type Index sh :: *
   -- Ix.range
   indices :: sh -> [Index sh]
   -- Ix.index
   offset :: sh -> Index sh -> Int
   offset sh = snd $ sizeOffset sh
   -- Ix.unsafeIndex
   uncheckedOffset :: sh -> Index sh -> Int
   uncheckedOffset = offset
   -- Ix.inRange
   inBounds :: sh -> Index sh -> Bool

   sizeOffset :: sh -> (Int, Index sh -> Int)
   sizeOffset sh = (size sh, offset sh)
   uncheckedSizeOffset :: sh -> (Int, Index sh -> Int)
   uncheckedSizeOffset sh = (uncheckedSize sh, uncheckedOffset sh)

class Indexed sh => InvIndexed sh where
   {- |
   It should hold @indexFromOffset sh k == indices sh !! k@,
   but 'indexFromOffset' should generally be faster.
   -}
   indexFromOffset :: sh -> Int -> Index sh
   uncheckedIndexFromOffset :: sh -> Int -> Index sh
   uncheckedIndexFromOffset = indexFromOffset

errorIndexFromOffset :: String -> Int -> a
errorIndexFromOffset name k =
   error $ printf "indexFromOffset (%s): index %d out of range" name k


instance C () where
   size () = 1
   uncheckedSize () = 1

instance Indexed () where
   type Index () = ()
   indices () = [()]
   offset () () = 0
   uncheckedOffset () () = 0
   inBounds () () = True

instance InvIndexed () where
   indexFromOffset () 0 = ()
   indexFromOffset () k = errorIndexFromOffset "()" k
   uncheckedIndexFromOffset () _ = ()


{- |
'ZeroBased' denotes a range starting at zero and has a certain length.
-}
newtype ZeroBased n = ZeroBased {zeroBasedSize :: n}
   deriving (Eq, Show)

instance Functor ZeroBased where
   fmap f (ZeroBased n) = ZeroBased $ f n

instance Applicative ZeroBased where
   pure = ZeroBased
   ZeroBased f <*> ZeroBased n = ZeroBased $ f n

instance (NFData n) => NFData (ZeroBased n) where
   rnf (ZeroBased n) = rnf n

instance (Storable n) => Storable (ZeroBased n) where
   sizeOf = Store.sizeOf zeroBasedSize
   alignment = Store.alignment zeroBasedSize
   peek = Store.peek ZeroBased
   poke = Store.poke zeroBasedSize

instance (Integral n) => C (ZeroBased n) where
   size (ZeroBased len) = fromIntegral len
   uncheckedSize (ZeroBased len) = fromIntegral len

instance (Integral n) => Indexed (ZeroBased n) where
   type Index (ZeroBased n) = n
   indices (ZeroBased len) = indices $ Shifted 0 len
   offset (ZeroBased len) = offset $ Shifted 0 len
   uncheckedOffset _ ix = fromIntegral ix
   inBounds (ZeroBased len) ix = 0<=ix && ix<len

instance (Integral n) => InvIndexed (ZeroBased n) where
   indexFromOffset (ZeroBased len) k0 =
      let k = fromIntegral k0
      in  if 0<=k && k<len
            then k
            else errorIndexFromOffset "ZeroBased" k0
   uncheckedIndexFromOffset _ k = fromIntegral k


{- |
'OneBased' denotes a range starting at one and has a certain length.
-}
newtype OneBased n = OneBased {oneBasedSize :: n}
   deriving (Eq, Show)

instance Functor OneBased where
   fmap f (OneBased n) = OneBased $ f n

instance Applicative OneBased where
   pure = OneBased
   OneBased f <*> OneBased n = OneBased $ f n

instance (NFData n) => NFData (OneBased n) where
   rnf (OneBased n) = rnf n

instance (Storable n) => Storable (OneBased n) where
   sizeOf = Store.sizeOf oneBasedSize
   alignment = Store.alignment oneBasedSize
   peek = Store.peek OneBased
   poke = Store.poke oneBasedSize

instance (Integral n) => C (OneBased n) where
   size (OneBased len) = fromIntegral len
   uncheckedSize (OneBased len) = fromIntegral len

instance (Integral n) => Indexed (OneBased n) where
   type Index (OneBased n) = n
   indices (OneBased len) = indices $ Shifted 1 len
   offset (OneBased len) = offset $ Shifted 1 len
   uncheckedOffset _ ix = fromIntegral ix - 1
   inBounds (OneBased len) ix = 0<ix && ix<=len

instance (Integral n) => InvIndexed (OneBased n) where
   indexFromOffset (OneBased len) k0 =
      let k = fromIntegral k0
      in  if 0<=k && k<len
            then 1+k
            else errorIndexFromOffset "OneBased" k0
   uncheckedIndexFromOffset _ k = 1 + fromIntegral k


{- |
'Range' denotes an inclusive range like
those of the Haskell 98 standard @Array@ type from the @array@ package.
E.g. the shape type @(Range Int32, Range Int64)@
is equivalent to the ix type @(Int32, Int64)@ for @Array@s.
-}
data Range n = Range {rangeFrom, rangeTo :: n}
   deriving (Eq, Show)

instance Functor Range where
   fmap f (Range from to) = Range (f from) (f to)

instance (NFData n) => NFData (Range n) where
   rnf (Range from to) = rnf (from,to)

instance (Ix.Ix n) => C (Range n) where
   size (Range from to) = Ix.rangeSize (from,to)
   uncheckedSize (Range from to) = Ix.unsafeRangeSize (from,to)

instance (Ix.Ix n) => Indexed (Range n) where
   type Index (Range n) = n
   indices (Range from to) = Ix.range (from,to)
   offset (Range from to) ix = Ix.index (from,to) ix
   uncheckedOffset (Range from to) ix = Ix.unsafeIndex (from,to) ix
   inBounds (Range from to) ix = Ix.inRange (from,to) ix

-- pretty inefficient when we rely solely on Ix
instance (Ix.Ix n) => InvIndexed (Range n) where
   indexFromOffset (Range from to) k =
      if 0<=k && k < Ix.rangeSize (from,to)
         then Ix.range (from,to) !! k
         else errorIndexFromOffset "Range" k
   uncheckedIndexFromOffset (Range from to) k = Ix.range (from,to) !! k

-- cf. sample-frame:Stereo
instance Storable n => Storable (Range n) where
   {-# INLINE sizeOf #-}
   {-# INLINE alignment #-}
   {-# INLINE peek #-}
   {-# INLINE poke #-}
   sizeOf ~(Range l r) = sizeOf l + mod (- sizeOf l) (alignment r) + sizeOf r
   alignment ~(Range l _) = alignment l
   poke p (Range l r) =
      let q = castToElemPtr p
      in  poke q l >> pokeElemOff q 1 r
   peek p =
      let q = castToElemPtr p
      in  Monad.lift2 Range (peek q) (peekElemOff q 1)


{- |
'Shifted' denotes a range defined by the start index and the length.
-}
data Shifted n = Shifted {shiftedOffset, shiftedSize :: n}
   deriving (Eq, Show)

instance Functor Shifted where
   fmap f (Shifted from to) = Shifted (f from) (f to)

instance (NFData n) => NFData (Shifted n) where
   rnf (Shifted from to) = rnf (from,to)

instance (Integral n) => C (Shifted n) where
   size (Shifted _offs len) = fromIntegral len
   uncheckedSize (Shifted _offs len) = fromIntegral len

instance (Integral n) => Indexed (Shifted n) where
   type Index (Shifted n) = n
   indices (Shifted offs len) =
      map snd $
      takeWhile ((>0) . fst) $
      zip
         (iterate (subtract 1) len)
         (iterate (1+) offs)
   offset (Shifted offs len) ix =
      if ix<offs
        then error "Shape.Shifted: array index too small"
        else
          let k = ix-offs
          in  if k<len
                then fromIntegral k
                else error "Shape.Shifted: array index too big"
   uncheckedOffset (Shifted offs _len) ix = fromIntegral $ ix-offs
   inBounds (Shifted offs len) ix = ix < offs+len

instance (Integral n) => InvIndexed (Shifted n) where
   indexFromOffset (Shifted offs len) k0 =
      let k = fromIntegral k0
      in  if 0<=k && k<len
            then offs+k
            else errorIndexFromOffset "Shifted" k0
   uncheckedIndexFromOffset (Shifted offs _len) k = offs + fromIntegral k

-- cf. sample-frame:Stereo
instance Storable n => Storable (Shifted n) where
   {-# INLINE sizeOf #-}
   {-# INLINE alignment #-}
   {-# INLINE peek #-}
   {-# INLINE poke #-}
   sizeOf ~(Shifted l n) = sizeOf l + mod (- sizeOf l) (alignment n) + sizeOf n
   alignment ~(Shifted l _) = alignment l
   poke p (Shifted l n) =
      let q = castToElemPtr p
      in  poke q l >> pokeElemOff q 1 n
   peek p =
      let q = castToElemPtr p
      in  Monad.lift2 Shifted (peek q) (peekElemOff q 1)


{-# INLINE castToElemPtr #-}
castToElemPtr :: Ptr (f a) -> Ptr a
castToElemPtr = castPtr



{- |
'Enumeration' denotes a shape of fixed size
that is defined by 'Enum' and 'Bounded' methods.
For correctness it is necessary that the 'Enum' and 'Bounded'
are properly implemented.
Automatically derived instances are fine.
-}
data Enumeration n = Enumeration
   deriving (Eq, Show)

instance NFData (Enumeration n) where
   rnf Enumeration = ()

instance (Enum n, Bounded n) => C (Enumeration n) where
   size = uncheckedSize
   uncheckedSize sh = intFromEnum sh maxBound - intFromEnum sh minBound + 1

instance (Enum n, Bounded n) => Indexed (Enumeration n) where
   type Index (Enumeration n) = n
   indices sh = [asEnumType sh minBound .. asEnumType sh maxBound]
   offset = uncheckedOffset
   uncheckedOffset sh ix = fromEnum ix - intFromEnum sh minBound
   inBounds _sh _ix = True

instance (Enum n, Bounded n) => InvIndexed (Enumeration n) where
   indexFromOffset sh k =
      if 0<=k && k <= intFromEnum sh maxBound - intFromEnum sh minBound
         then uncheckedIndexFromOffset sh k
         else errorIndexFromOffset "Enumeration" k
   uncheckedIndexFromOffset sh k = toEnum $ intFromEnum sh minBound + k

asEnumType :: Enumeration n -> n -> n
asEnumType Enumeration = id

intFromEnum :: (Enum n) => Enumeration n -> n -> Int
intFromEnum Enumeration = fromEnum

instance Storable (Enumeration n) where
   {-# INLINE sizeOf #-}
   {-# INLINE alignment #-}
   {-# INLINE peek #-}
   {-# INLINE poke #-}
   sizeOf ~Enumeration = 0
   alignment ~Enumeration = 1
   poke _p Enumeration = return ()
   peek _p = return Enumeration


{- |
This data type wraps another array shape.
Its index type is a wrapped 'Int'.
The advantages are:
No conversion forth and back 'Int' and @Index sh@.
You can convert once using 'deferIndex' and 'revealIndex'
whenever you need your application specific index type.
No need for e.g. @Storable (Index sh)@, because 'Int' is already 'Storable'.
-}
newtype Deferred sh = Deferred sh
   deriving (Eq, Show)

newtype DeferredIndex ix = DeferredIndex Int
   deriving (Eq, Show)

instance (NFData sh) => NFData (Deferred sh) where
   rnf (Deferred sh) = rnf sh

instance (C sh) => C (Deferred sh) where
   size (Deferred sh) = size sh
   uncheckedSize (Deferred sh) = uncheckedSize sh

instance (C sh) => Indexed (Deferred sh) where
   type Index (Deferred sh) = DeferredIndex (Index sh)
   indices (Deferred sh) = map DeferredIndex $ take (size sh) [0 ..]
   offset (Deferred sh) (DeferredIndex k) = offset (ZeroBased $ size sh) k
   uncheckedOffset _ (DeferredIndex k) = k
   inBounds (Deferred sh) (DeferredIndex k) =
      inBounds (ZeroBased $ size sh) k

instance (C sh) => InvIndexed (Deferred sh) where
   indexFromOffset sh k =
      DeferredIndex $ indexFromOffset (ZeroBased $ size sh) k
   uncheckedIndexFromOffset _sh = DeferredIndex

deferIndex :: (Indexed sh, Index sh ~ ix) => sh -> ix -> DeferredIndex ix
deferIndex sh ix = DeferredIndex $ offset sh ix

revealIndex :: (InvIndexed sh, Index sh ~ ix) => sh -> DeferredIndex ix -> ix
revealIndex sh (DeferredIndex ix) = indexFromOffset sh ix

instance Storable (DeferredIndex ix) where
   {-# INLINE sizeOf #-}
   {-# INLINE alignment #-}
   {-# INLINE peek #-}
   {-# INLINE poke #-}
   sizeOf (DeferredIndex k) = sizeOf k
   alignment (DeferredIndex k) = alignment k
   poke p (DeferredIndex k) = poke (castPtr p) k
   peek p = fmap DeferredIndex $ peek (castPtr p)



{- |
Row-major composition of two dimensions.
-}
instance (C sh0, C sh1) => C (sh0,sh1) where
   size (sh0,sh1) = size sh0 * size sh1
   uncheckedSize (sh0,sh1) = uncheckedSize sh0 * uncheckedSize sh1

instance (Indexed sh0, Indexed sh1) => Indexed (sh0,sh1) where
   type Index (sh0,sh1) = (Index sh0, Index sh1)
   indices (sh0,sh1) = Monad.lift2 (,) (indices sh0) (indices sh1)
   offset (sh0,sh1) =
      offset sh0 . fst
      `combineOffset`
      mapSnd (.snd) (sizeOffset sh1)
   uncheckedOffset (sh0,sh1) =
      uncheckedOffset sh0 . fst
      `combineOffset`
      mapSnd (.snd) (uncheckedSizeOffset sh1)
   sizeOffset (sh0,sh1) =
      mapSnd (.fst) (sizeOffset sh0)
      `combineSizeOffset`
      mapSnd (.snd) (sizeOffset sh1)
   uncheckedSizeOffset (sh0,sh1) =
      mapSnd (.fst) (uncheckedSizeOffset sh0)
      `combineSizeOffset`
      mapSnd (.snd) (uncheckedSizeOffset sh1)
   inBounds (sh0,sh1) (ix0,ix1) = inBounds sh0 ix0 && inBounds sh1 ix1

instance (InvIndexed sh0, InvIndexed sh1) => InvIndexed (sh0,sh1) where
   indexFromOffset (sh0,sh1) k =
      runInvIndex k $ liftA2 (,) (pickLastIndex sh0) (pickIndex sh1)
   uncheckedIndexFromOffset (sh0,sh1) k =
      runInvIndex k $ liftA2 (,) (uncheckedPickLastIndex sh0) (pickIndex sh1)


instance (C sh0, C sh1, C sh2) => C (sh0,sh1,sh2) where
   size (sh0,sh1,sh2) = size sh0 * size sh1 * size sh2
   uncheckedSize (sh0,sh1,sh2) =
      uncheckedSize sh0 * uncheckedSize sh1 * uncheckedSize sh2

instance (Indexed sh0, Indexed sh1, Indexed sh2) => Indexed (sh0,sh1,sh2) where
   type Index (sh0,sh1,sh2) = (Index sh0, Index sh1, Index sh2)
   indices (sh0,sh1,sh2) =
      Monad.lift3 (,,) (indices sh0) (indices sh1) (indices sh2)
   uncheckedOffset (sh0,sh1,sh2) =
      uncheckedOffset sh0 . fst3
      `combineOffset`
      mapSnd (.snd3) (uncheckedSizeOffset sh1)
      `combineSizeOffset`
      mapSnd (.thd3) (uncheckedSizeOffset sh2)
   sizeOffset (sh0,sh1,sh2) =
      mapSnd (.fst3) (sizeOffset sh0)
      `combineSizeOffset`
      mapSnd (.snd3) (sizeOffset sh1)
      `combineSizeOffset`
      mapSnd (.thd3) (sizeOffset sh2)
   uncheckedSizeOffset (sh0,sh1,sh2) =
      mapSnd (.fst3) (uncheckedSizeOffset sh0)
      `combineSizeOffset`
      mapSnd (.snd3) (uncheckedSizeOffset sh1)
      `combineSizeOffset`
      mapSnd (.thd3) (uncheckedSizeOffset sh2)
   inBounds (sh0,sh1,sh2) (ix0,ix1,ix2) =
      inBounds sh0 ix0 && inBounds sh1 ix1 && inBounds sh2 ix2

instance
   (InvIndexed sh0, InvIndexed sh1, InvIndexed sh2) =>
      InvIndexed (sh0,sh1,sh2) where
   indexFromOffset (sh0,sh1,sh2) k =
      runInvIndex k $
      liftA3 (,,) (pickLastIndex sh0) (pickIndex sh1) (pickIndex sh2)
   uncheckedIndexFromOffset (sh0,sh1,sh2) k =
      runInvIndex k $
      liftA3 (,,) (uncheckedPickLastIndex sh0) (pickIndex sh1) (pickIndex sh2)

runInvIndex :: s -> Back.Backwards (MS.State s) a -> a
runInvIndex k = flip MS.evalState k . Back.forwards

pickLastIndex ::
   (InvIndexed sh) => sh -> Back.Backwards (MS.State Int) (Index sh)
pickLastIndex sh =
   Back.Backwards $ MS.gets $ indexFromOffset sh

uncheckedPickLastIndex ::
   (InvIndexed sh) => sh -> Back.Backwards (MS.State Int) (Index sh)
uncheckedPickLastIndex sh =
   Back.Backwards $ MS.gets $ uncheckedIndexFromOffset sh

pickIndex :: (InvIndexed sh) => sh -> Back.Backwards (MS.State Int) (Index sh)
pickIndex sh =
   fmap (uncheckedIndexFromOffset sh) $
   Back.Backwards $ MS.state $ \k -> swap $ divMod k $ size sh



infixr 7 `combineOffset`, `combineSizeOffset`

{-# INLINE combineOffset #-}
combineOffset :: Num a => (ix -> a) -> (a, ix -> a) -> ix -> a
combineOffset offset0 (size1,offset1) ix = offset0 ix * size1 + offset1 ix

{-# INLINE combineSizeOffset #-}
combineSizeOffset :: Num a => (a, ix -> a) -> (a, ix -> a) -> (a, ix -> a)
combineSizeOffset (size0,offset0) (size1,offset1) =
   (size0*size1, \ix -> offset0 ix * size1 + offset1 ix)



data Lower = Lower deriving (Eq, Show)
data Upper = Upper deriving (Eq, Show)

class TriangularPart part where
   switchTriangularPart :: f Lower -> f Upper -> f part
instance TriangularPart Lower where switchTriangularPart f _ = f
instance TriangularPart Upper where switchTriangularPart _ f = f

getConstAs :: c -> Const a c -> a
getConstAs _ = getConst

caseTriangularPart :: (TriangularPart part) => part -> a -> a -> a
caseTriangularPart part lo up =
   getConstAs part $ switchTriangularPart (Const lo) (Const up)

data Triangular part size =
   Triangular {
      triangularPart :: part,
      triangularSize :: size
   } deriving (Eq, Show)

type LowerTriangular = Triangular Lower
type UpperTriangular = Triangular Upper

lowerTriangular :: size -> LowerTriangular size
lowerTriangular = Triangular Lower

upperTriangular :: size -> UpperTriangular size
upperTriangular = Triangular Upper

-- cf. Data.Bifunctor.Flip
newtype Flip f b a = Flip {getFlip :: f a b}

instance
      (TriangularPart part, NFData size) => NFData (Triangular part size) where
   rnf (Triangular part sz) =
      rnf
         (flip getFlip part $
            switchTriangularPart (Flip $ \Lower -> ()) (Flip $ \Upper -> ()),
          sz)

instance (TriangularPart part, C size) => C (Triangular part size) where
   size (Triangular _part sz) = triangleSize $ size sz
   uncheckedSize (Triangular _part sz) = triangleSize $ uncheckedSize sz

instance
   (TriangularPart part, Indexed size) =>
      Indexed (Triangular part size) where
   type Index (Triangular part size) = (Index size, Index size)

   indices (Triangular part sz) =
      let ixs = indices sz
      in concat $
         caseTriangularPart part
            (zipWith (\cs r -> map ((,) r) cs)
               (NonEmpty.tail $ NonEmpty.inits ixs) ixs)
            (zipWith (\r cs -> map ((,) r) cs) ixs $ tails ixs)

   uncheckedOffset sh = snd $ uncheckedSizeOffset sh

   sizeOffset (Triangular part sz) =
      let (n, getOffset) = sizeOffset sz
      in (triangleSize n, \(rs,cs) ->
            let r = getOffset rs
                c = getOffset cs
            in if compareIndices part r c
                  then triangleOffset part n (r,c)
                  else error "Shape.Triangular.sizeOffset: wrong array part")

   uncheckedSizeOffset (Triangular part sz) =
      let (n, getOffset) = uncheckedSizeOffset sz
      in (triangleSize n, \(rs,cs) ->
            triangleOffset part n (getOffset rs, getOffset cs))

   inBounds (Triangular part sz) ix@(r,c) =
      inBounds (sz,sz) ix
      &&
      let getOffset = offset sz
      in compareIndices part (getOffset r) (getOffset c)

triangleOffset :: TriangularPart part => part -> Int -> (Int, Int) -> Int
triangleOffset part n (r,c) =
   caseTriangularPart part
      (triangleSize r + c)
      (triangleSize n - triangleSize (n-r) + c-r)

compareIndices :: (TriangularPart part, Ord a) => part -> a -> a -> Bool
compareIndices part = caseTriangularPart part (>=) (<=)

instance
   (TriangularPart part, InvIndexed size) =>
      InvIndexed (Triangular part size) where

   indexFromOffset (Triangular part sz) k =
      mapPair (indexFromOffset sz, indexFromOffset sz) $
      caseTriangularPart part
         (let r = floor (triangleRootDouble k)
          in (r, k - triangleSize r))
         (let n = size sz
              triSize = triangleSize n
              rr = ceiling (triangleRootDouble (triSize-k))
              r = n - rr
          in (r, k+r - (triSize - triangleSize rr)))

triangleSize :: Int -> Int
triangleSize n = div (n*(n+1)) 2

{-
n*(n+1)/2 = m
n^2 + n - 2m = 0
n = -1/2 + sqrt(1/4+2m)
  = (sqrt(8m+1) - 1) / 2
-}
triangleRoot :: Floating a => a -> a
triangleRoot sz = (sqrt (8*sz+1)-1)/2

triangleRootDouble :: Int -> Double
triangleRootDouble = triangleRoot . fromIntegral



infixr 5 :+:

data sh0:+:sh1 = sh0:+:sh1
   deriving (Eq, Show)

instance (NFData sh0, NFData sh1) => NFData (sh0:+:sh1) where
   rnf (sh0:+:sh1) = rnf (sh0,sh1)

instance (C sh0, C sh1) => C (sh0:+:sh1) where
   size (sh0:+:sh1) = size sh0 + size sh1
   uncheckedSize (sh0:+:sh1) = uncheckedSize sh0 + uncheckedSize sh1

instance (Indexed sh0, Indexed sh1) => Indexed (sh0:+:sh1) where
   type Index (sh0:+:sh1) = Either (Index sh0) (Index sh1)
   indices (sh0:+:sh1) = map Left (indices sh0) ++ map Right (indices sh1)
   offset (sh0:+:sh1) ix =
      case ix of
         Left ix0 -> offset sh0 ix0
         Right ix1 -> size sh0 + offset sh1 ix1
   uncheckedOffset (sh0:+:sh1) ix =
      case ix of
         Left ix0 -> uncheckedOffset sh0 ix0
         Right ix1 -> uncheckedSize sh0 + uncheckedOffset sh1 ix1
   sizeOffset (sh0:+:sh1) =
      let (n0, getOffset0) = sizeOffset sh0
          (n1, getOffset1) = sizeOffset sh1
      in (n0+n1, either getOffset0 ((n0+) . getOffset1))
   uncheckedSizeOffset (sh0:+:sh1) =
      let (n0, getOffset0) = uncheckedSizeOffset sh0
          (n1, getOffset1) = uncheckedSizeOffset sh1
      in (n0+n1, either getOffset0 ((n0+) . getOffset1))
   inBounds (sh0:+:sh1) = either (inBounds sh0) (inBounds sh1)

instance (InvIndexed sh0, InvIndexed sh1) => InvIndexed (sh0:+:sh1) where
   indexFromOffset (sh0:+:sh1) k =
      let pivot = size sh0
      in if k < pivot
            then Left $ indexFromOffset sh0 k
            else Right $ indexFromOffset sh1 $ k-pivot
   uncheckedIndexFromOffset (sh0:+:sh1) k =
      let pivot = size sh0
      in if k < pivot
            then Left $ uncheckedIndexFromOffset sh0 k
            else Right $ uncheckedIndexFromOffset sh1 $ k-pivot
