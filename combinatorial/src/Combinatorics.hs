{- |
Count and create combinatorial objects.
Also see 'combinat' package.
-}
module Combinatorics (
   permute,
   permuteFast,
   permuteShare,
   permuteRep,
   choose,
   variateRep,
   variate,
   tuples,
   partitions,
   rectifications,
   setPartitions,
   chooseUnrank,
   chooseUnrankMaybe,
   chooseRank,
   factorial,
   binomial,
   binomialSeq,
   binomialGen,
   binomialSeqGen,
   multinomial,
   factorials,
   binomials,
   catalanNumber,
   catalanNumbers,
   derangementNumber,
   derangementNumbers,
   setPartitionNumbers,
   surjectiveMappingNumber,
   surjectiveMappingNumbers,
   fibonacciNumber,
   fibonacciNumbers,
   ) where

import qualified PowerSeries
import qualified Combinatorics.Private as Comb

import Data.Function.HT (nest, )
import Data.Maybe.HT (toMaybe, )
import Data.Tuple.HT (mapFst, )
import qualified Data.List.Match as Match
import Data.List.HT (mapAdjacent, removeEach, )
import Data.List (genericIndex, )

import Control.Monad (liftM2, )


{-* Generate compositions from a list of elements. -}

-- several functions for permutation
-- cf. Equation.hs

{- |
Generate list of all permutations of the input list.
The list is sorted lexicographically.
-}
permute :: [a] -> [[a]]
permute = Comb.permuteRec

{- |
Generate list of all permutations of the input list.
It is not lexicographically sorted.
It is slightly faster and consumes less memory
than the lexicographical ordering 'permute'.
-}
permuteFast :: [a] -> [[a]]
permuteFast x = permuteFastStep [] x []

{- |
Each element of (allcycles x) has a different element at the front.
Iterate cycling on the tail elements of each element list of (allcycles x).
-}
permuteFastStep :: [a] -> [a] -> [[a]] -> [[a]]
permuteFastStep suffix [] tl = suffix:tl
permuteFastStep suffix x  tl =
   foldr (\c -> permuteFastStep (head c : suffix) (tail c)) tl (allCycles x)

{- |
All permutations share as much suffixes as possible.
The reversed permutations are sorted lexicographically.
-}
permuteShare :: [a] -> [[a]]
permuteShare x =
   map fst $
--   map (\(y,[]) -> y) $  -- safer but inefficient
   nest (length x) (concatMap permuteShareStep) [([], x)]

permuteShareStep :: ([a], [a]) -> [([a], [a])]
permuteShareStep (perm,todo) =
   map
      (mapFst (:perm))
      (removeEach todo)


permuteRep :: [(a,Int)] -> [[a]]
permuteRep = Comb.permuteRep


choose :: Int -> Int -> [[Bool]]
choose = Comb.chooseRec


{- |
Generate all choices of n elements out of the list x with repetitions.
\"variation\" seems to be used historically,
but I like it more than \"k-permutation\".
-}
variateRep :: Int -> [a] -> [[a]]
variateRep = Comb.variateRep


{- |
Generate all choices of n elements out of the list x without repetitions.
It holds
   @ variate (length xs) xs == permute xs @
-}
variate :: Int -> [a] -> [[a]]
variate = Comb.variateRec


{- |
Generate all choices of n elements out of the list x
respecting the order in x and without repetitions.
-}
tuples :: Int -> [a] -> [[a]]
tuples = Comb.tuplesRec


partitions :: [a] -> [([a],[a])]
partitions =
   foldr
      (\x -> concatMap (\(lxs,rxs) -> [(x:lxs,rxs), (lxs,x:rxs)]))
      [([],[])]

{- |
Number of possibilities arising in rectification of a predicate
in deductive database theory.
Stefan Brass, \"Logische Programmierung und deduktive Datenbanken\", 2007,
page 7-60
This is isomorphic to the partition of @n@-element sets
into @k@ non-empty subsets.
<http://oeis.org/A048993>

> *Combinatorics> map (length . uncurry rectifications) $ do x<-[0..10]; y<-[0..x]; return (x,[1..y::Int])
> [1,0,1,0,1,1,0,1,3,1,0,1,7,6,1,0,1,15,25,10,1,0,1,31,90,65,15,1,0,1,63,301,350,140,21,1,0,1,127,966,1701,1050,266,28,1,0,1,255,3025,7770,6951,2646,462,36,1,0,1,511,9330,34105,42525,22827,5880,750,45,1]
-}
rectifications :: Int -> [a] -> [[a]]
rectifications =
   let recourse _ 0 xt =
          if null xt
            then [[]]
            else []
       recourse ys n xt =
          let n1 = pred n
          in  liftM2 (:) ys (recourse ys n1 xt) ++
              case xt of
                 [] -> []
                 (x:xs) -> map (x:) (recourse (ys++[x]) n1 xs)
   in  recourse []

{- |
Their number is @k^n@.
-}
{-
setPartitionsEmpty :: Int -> [a] -> [[[a]]]
setPartitionsEmpty k =
   let recourse [] = [replicate k []]
       recourse (x:xs) =
          map (\(ys0,y,ys1) -> ys0 ++ [x:y] ++ ys1) $
          concatMap splitEverywhere (recourse xs)
{-
          do xs1 <- recourse xs
             (ys0,y,ys1) <- splitEverywhere xs1
             return (ys0 ++ [x:y] ++ ys1)
-}
   in  recourse
-}

setPartitions :: Int -> [a] -> [[[a]]]
setPartitions 0 xs =
   if null xs
     then [[]]
     else [  ]
setPartitions _ [] = []
setPartitions 1 xs = [[xs]]  -- unnecessary for correctness, but useful for efficiency
setPartitions k (x:xs) =
   do (rest, choosen) <- partitions xs
      part <- setPartitions (pred k) rest
      return ((x:choosen) : part)


{-* Rank and unrank combinatorial objects. -}

{- |
@chooseUnrank n k i == choose n k !! i@
-}
chooseUnrank :: Integral a => a -> a -> a -> [Bool]
chooseUnrank = Comb.chooseUnrankRec

chooseUnrankMaybe :: Int -> Int -> Int -> Maybe [Bool]
chooseUnrankMaybe n k i =
   toMaybe
      (0 <= i && i < binomial n k)
      (chooseUnrank n k i)
-- error ("chooseUnrank: out of range " ++ show (n, k, i))


{- |
<https://en.wikipedia.org/wiki/Combinatorial_number_system>
-}
chooseRank :: Integral a => [Bool] -> (a, a, a)
chooseRank =
   foldl
      (\(n,k0,i0) (bins,b) ->
        let (k1,i1) = if b then (succ k0, i0 + genericIndex (bins++[0]) k1) else (k0,i0)
        in  (succ n, k1, i1))
      (0,0,0) .
   zip binomials .
   reverse


{-* Generate complete lists of combinatorial numbers. -}

factorial :: Integral a => a -> a
factorial n = product [1..n]

{-| Pascal's triangle containing the binomial coefficients. -}
binomial :: Integral a => a -> a -> a
binomial = Comb.binomial

binomialSeq :: Integral a => a -> [a]
binomialSeq = Comb.binomialSeq


binomialGen :: (Integral a, Fractional b) => b -> a -> b
binomialGen n k = genericIndex (binomialSeqGen n) k

binomialSeqGen :: (Fractional b) => b -> [b]
binomialSeqGen n =
   scanl (\acc (num,den) -> acc*num / den) 1
         (zip (iterate (subtract 1) n) (iterate (1+) 1))


multinomial :: Integral a => [a] -> a
multinomial =
   product . mapAdjacent binomial . scanr1 (+)


{-* Generate complete lists of factorial numbers. -}

factorials :: Num a => [a]
factorials = Comb.factorials

{-|
Pascal's triangle containing the binomial coefficients.
Only efficient if a prefix of all rows is required.
It is not efficient for picking particular rows
or even particular elements.
-}
binomials :: Num a => [[a]]
binomials = Comb.binomials


{- |
@catalanNumber n@ computes the number of binary trees with @n@ nodes.
-}
catalanNumber :: Integer -> Integer
catalanNumber n =
   let (c,r) = divMod (binomial (2*n) n) (n+1)
   in  if r==0
         then c
         else error "catalanNumber: Integer implementation broken"

{- |
Compute the sequence of Catalan numbers by recurrence identity.
It is @catalanNumbers !! n == catalanNumber n@
-}
catalanNumbers :: Num a => [a]
catalanNumbers =
   let xs = 1 : PowerSeries.mul xs xs
   in  xs



derangementNumber :: Integer -> Integer
derangementNumber n =
   sum (scanl (*) ((-1) ^ mod n 2) [-n,1-n..(-1)])

{- |
Number of fix-point-free permutations with @n@ elements.

<http://oeis.org/A000166>
-}
derangementNumbers :: Num a => [a]
derangementNumbers = Comb.derangementNumbersPS0


-- generation of all possibilities and computation of their number should be in different modules

{- |
Number of partitions of an @n@ element set into @k@ non-empty subsets.
Known as Stirling numbers <http://oeis.org/A048993>.
-}
setPartitionNumbers :: Num a => [[a]]
setPartitionNumbers = Comb.setPartitionNumbers


{- |
@surjectiveMappingNumber n k@ computes the number of surjective mappings
from a @n@ element set to a @k@ element set.

<http://oeis.org/A019538>
-}
surjectiveMappingNumber :: Integer -> Integer -> Integer
surjectiveMappingNumber n k =
   foldl subtract 0 $
   zipWith (*)
      (map (^n) [0..])
      (binomialSeq k)

surjectiveMappingNumbers :: Num a => [[a]]
surjectiveMappingNumbers = Comb.surjectiveMappingNumbersPS


{- |
Multiply two Fibonacci matrices, that is matrices of the form

> /F[n-1] F[n]  \
> \F[n]   F[n+1]/
-}
fiboMul ::
   (Integer,Integer,Integer) ->
   (Integer,Integer,Integer) ->
   (Integer,Integer,Integer)
fiboMul (f0,f1,f2) (g0,g1,g2) =
   let h0 = f0*g0 + f1*g1
       h1 = f0*g1 + f1*g2
--     h1 = f1*g0 + f2*g1
       h2 = f1*g1 + f2*g2
   in  (h0,h1,h2)


{-
Fast computation using matrix power of

> /0 1\
> \1 1/

Hard-coded fast power with integer exponent.
Better use a generic algorithm.
-}
fibonacciNumber :: Integer -> Integer
fibonacciNumber x =
   let aux   0  = (1,0,1)
       aux (-1) = (-1,1,0)
       aux n =
          let (m,r) = divMod n 2
              f = aux m
              f2 = fiboMul f f
          in  if r==0
                then f2
                else fiboMul (0,1,1) f2
       (_,y,_) = aux x
   in  y


{- |
Number of possibilities to compose a 2 x n rectangle of n bricks.

>  |||   |--   --|
>  |||   |--   --|
-}
fibonacciNumbers :: [Integer]
fibonacciNumbers =
   let xs = 0 : ys
       ys = 1 : zipWith (+) xs ys
   in  xs



{- * Auxiliary functions -}

{- candidates for Useful -}

{- | Create a list of all possible rotations of the input list. -}
allCycles :: [a] -> [[a]]
allCycles x =
   Match.take x (map (Match.take x) (iterate tail (cycle x)))
