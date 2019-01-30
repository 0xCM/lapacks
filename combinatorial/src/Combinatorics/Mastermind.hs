module Combinatorics.Mastermind (
   Eval(..),
   evaluate,
   evaluateAll,
   formatEvalHistogram,
   numberDistinct,
   ) where

import qualified Combinatorics.Permutation.WithoutSomeFixpoints as PermWOFP
import Combinatorics (binomial)

import Text.Printf (printf)

import qualified Data.Map as Map; import Data.Map (Map)
import qualified Data.Foldable as Fold
import qualified Data.List.HT as ListHT
import Data.Tuple.HT (mapPair)


{- |
Cf. @board-games@ package.
-}
data Eval = Eval {black, white :: Int}
   deriving (Eq, Ord, Show)

{- |
Given the code and a guess, compute the evaluation.
-}
evaluate :: (Ord a) => [a] -> [a] -> Eval
evaluate code attempt =
   uncurry Eval $
   mapPair
      (length,
       Fold.sum . uncurry (Map.intersectionWith min) .
       mapPair (histogram,histogram) . unzip) $
   ListHT.partition (uncurry (==)) $
   zip code attempt

{-
*Combinatorics.Mastermind> filter ((Eval 2 0 ==) . evaluate "aabbb") $ replicateM 5 ['a'..'c']
["aaaaa","aaaac","aaaca","aaacc","aacaa","aacac","aacca","aaccc","acbcc","accbc","acccb","cabcc","cacbc","caccb","ccbbc","ccbcb","cccbb"]
-}

evaluateAll :: (Ord a) => [[a]] -> [a] -> Map Eval Int
evaluateAll codes attempt = histogram $ map (evaluate attempt) codes

formatEvalHistogram :: Map Eval Int -> String
formatEvalHistogram m =
   let n = maximum $ map (\(Eval b w) -> b+w) $ Map.keys m
   in  unlines $
       zipWith
          (\b ->
             unwords .
             map (\w -> printf "%6d" $ Map.findWithDefault 0 (Eval b w) m))
          [0..] (reverse $ tail $ ListHT.inits [0..n])


histogram :: (Ord a) => [a] -> Map a Int
histogram  =  Map.fromListWith (+) . map (\a -> (a,1))


{- |
@numberDistinct n k b w@ computes the number of matching codes,
given that all codes have distinct symbols.
@n@ is the alphabet size, @k@ the width of the code,
@b@ the number of black evaluation sticks and
@w@ the number of white evaluation sticks.
-}
numberDistinct :: Int -> Int -> Int -> Int -> Integer
numberDistinct n k b w =
   binomial (toInteger k) (toInteger b)
   *
   numberDistinctWhite (n-b) (k-b) w

{- |
@numberDistinctWhite n k w == numberDistinct n k 0 w@
-}
numberDistinctWhite :: Int -> Int -> Int -> Integer
numberDistinctWhite n k w =
   let ni = toInteger n
       ki = toInteger k
       wi = toInteger w
   in  binomial ki wi * PermWOFP.numbers !! k !! w * binomial (ni-ki) (ki-wi)
