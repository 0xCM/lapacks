module Main (main) where

import qualified Combinatorics.Permutation.WithoutSomeFixpoints as PermWOFP
import qualified Combinatorics.Mastermind as Mastermind
import qualified Combinatorics.CardPairs as CardPairs
import qualified Combinatorics.Partitions as Parts
import qualified Combinatorics.BellNumbers as Bell
import qualified Combinatorics.Private as CombPriv
import qualified Combinatorics as Comb

import qualified Test.QuickCheck as QC
import Test.QuickCheck (Testable, quickCheck, )

import Control.Monad (liftM2, liftM3, )
import Control.Applicative ((<$>), )

import qualified Data.List.Match as Match
import qualified Data.List.Key as Key
import qualified Data.List as List
import qualified Data.Set as Set
import Data.Array ((!), )
import Data.Tuple.HT (uncurry3, )
import Data.List.HT (allEqual, isAscending, )
import Data.List (sort, nub, )
import Data.Eq.HT (equating, )



permuteSum :: [Int] -> Bool
permuteSum xs =
   sum (map sum (Comb.permute xs)) ==
   sum xs * Comb.factorial (length xs)

permute :: Ord a => [a] -> Bool
permute xs =
   allEqual $
   map (\p -> sort (p xs)) $
      Comb.permute :
      Comb.permuteFast :
      Comb.permuteShare :
      []

permuteAlt :: Eq a => [a] -> Bool
permuteAlt xs = CombPriv.permuteRec xs == CombPriv.permuteMSL xs


genPermuteRep :: Int -> QC.Gen [(Char, Int)]
genPermuteRep n = do
   xns <- QC.listOf $ liftM2 (,) QC.arbitrary $ QC.choose (0,n)
   return $ Match.take (takeWhile (<=n) $ scanl1 (+) $ map snd xns) xns

permuteRepM :: Eq a => [(a, Int)] -> Bool
permuteRepM xs = CombPriv.permuteRep xs == CombPriv.permuteRepM xs

permuteRepNub :: Eq a => [(a, Int)] -> Bool
permuteRepNub xs =
   let perms = Comb.permuteRep $ Key.nub fst xs
   in  perms == nub perms

permuteRepNubBig :: Ord a => [(a, Int)] -> Bool
permuteRepNubBig xs =
   let perms = Comb.permuteRep $ Key.nub fst xs
   in  List.sort perms == Set.toList (Set.fromList perms)

permuteRepMonotony :: Ord a => [(a, Int)] -> Bool
permuteRepMonotony = isAscending . Comb.permuteRep . Key.nub fst . sort

permuteRepChoose :: Int -> Int -> Bool
permuteRepChoose n k =
   Comb.choose n k == Comb.permuteRep [(False, n-k), (True, k)]

chooseLength :: Int -> Int -> Bool
chooseLength n k =
   all
      (\x  ->  n == length x  &&  k == length (filter id x))
      (Comb.choose n k)


genChoose :: QC.Gen (Int, Int)
genChoose = do
   n <- QC.choose (0,15)
   k <- QC.choose (-2,n)
   return (n,k)

choose :: Int -> Int -> Bool
choose n k =
   allEqual $
      CombPriv.chooseRec n k :
      CombPriv.chooseMSL n k :
      CombPriv.chooseMSL0 n k :
      []

genChooseIndex :: QC.Gen (Integer, Integer, Integer)
genChooseIndex = do
   n <- QC.choose (0,25)
   k <- QC.choose (0,n)
   i <- QC.choose (0, Comb.binomial n k - 1)
   return (n,k,i)

chooseUnrank :: Integer -> Integer -> Integer -> Bool
chooseUnrank n k i =
   CombPriv.chooseUnrankRec n k i  ==  CombPriv.chooseUnrankList n k i

chooseUnrankSequence :: Int -> Int -> Bool
chooseUnrankSequence n k =
   map (Comb.chooseUnrank n k) [0 .. Comb.binomial n k - 1]
     ==  Comb.choose n k

chooseRankUnrank :: Integer -> Integer -> Integer -> Bool
chooseRankUnrank n k i =
   Comb.chooseRank (Comb.chooseUnrank n k i)  ==  (n, k, i)

chooseUnrankRank :: [Bool] -> Bool
chooseUnrankRank bs =
   uncurry3 Comb.chooseUnrank
      (Comb.chooseRank bs :: (Integer, Integer, Integer))
     ==  bs



genVariate :: QC.Gen [Char]
genVariate = take 7 <$> QC.arbitrary

variate :: Eq a => Int -> [a] -> Bool
variate n xs =
   CombPriv.variateRec n xs == CombPriv.variateMSL n xs

variateRep :: Eq a => Int -> [a] -> Bool
variateRep n xs =
   CombPriv.variateRep n xs == CombPriv.variateRepM n xs

variatePermute :: Eq a => [a] -> Bool
variatePermute xs =
   Comb.variate (length xs) xs == Comb.permute xs

variatePermuteClip :: Eq a => Int -> [a] -> Bool
variatePermuteClip n xs =
   equating (take n) (Comb.variate (length xs) xs) (Comb.permute xs)



genTuples :: QC.Gen (Int, [Char])
genTuples = do
   xs <- take 16 <$> QC.arbitrary
   n <- QC.choose (-1, length xs + 1)
   return (n,xs)

tuples :: Eq a => Int -> [a] -> Bool
tuples n xs =
   allEqual $
      CombPriv.tuplesRec n xs :
      CombPriv.tuplesRec0 n xs :
      CombPriv.tuplesMSL n xs :
      CombPriv.tuplesMSL0 n xs :
      []


_setPartitionsMonotony :: Ord a => Int -> [a] -> Bool
_setPartitionsMonotony k =
   isAscending . Comb.setPartitions k . nub . sort

rectificationsMonotony :: Ord a => Int -> [a] -> Bool
rectificationsMonotony k =
   isAscending . Comb.rectifications k . nub . sort



factorial :: [Char] -> Bool
factorial xs =
   length (Comb.permute xs) == Comb.factorial (length xs)


binomial :: [Char] -> Int -> Bool
binomial xs k =
   length (Comb.tuples k xs) == Comb.binomial (length xs) k


genBinomial :: QC.Gen (Integer, Integer)
genBinomial = do
   n <- QC.choose (0,100)
   k <- QC.choose (0,n)
   return (n,k)

binomialFactorial :: Integer -> Integer -> Bool
binomialFactorial n k =
   let (q, r) =
         divMod
            (Comb.factorial n)
            (Comb.factorial k * Comb.factorial (n-k))
   in  r == 0 && Comb.binomial n k == q


binomialChoose :: Int -> Int -> Bool
binomialChoose n k =
   length (Comb.choose n k) == Comb.binomial n k

multinomialPermuteRep :: [(Char,Int)] -> Bool
multinomialPermuteRep xs =
   length (Comb.permuteRep xs) == Comb.multinomial (map snd xs)

multinomialCommutative :: [Integer] -> Bool
multinomialCommutative xs =
   Comb.multinomial xs == Comb.multinomial (sort xs)

setPartitionNumbers :: Int -> [Int] -> Bool
setPartitionNumbers k xs =
   length (Comb.setPartitions k xs) ==
   (Comb.setPartitionNumbers !! length xs ++ repeat 0) !! k

rectificationNumbers :: Int -> [Int] -> Bool
rectificationNumbers k xs =
   length (Comb.rectifications k xs) ==
   (Comb.setPartitionNumbers !! k ++ repeat 0) !! length xs


surjectiveMappingNumber :: Int -> Bool
surjectiveMappingNumber =
   equalFuncList2 Comb.surjectiveMappingNumber Comb.surjectiveMappingNumbers

surjectiveMappingNumbers :: Int -> Bool
surjectiveMappingNumbers n =
   allEqual $ map (take n) $ (
      CombPriv.surjectiveMappingNumbersPS :
      CombPriv.surjectiveMappingNumbersStirling :
      [] :: [[[Integer]]])


equalFuncList :: (Integer -> Integer) -> [Integer] -> Int -> Bool
equalFuncList f xs n =
   equating (take n) xs (map f $ iterate (1+) 0)

factorials :: Int -> Bool
factorials = equalFuncList Comb.factorial Comb.factorials

equalFuncList2 :: (Integer -> Integer -> Integer) -> [[Integer]] -> Int -> Bool
equalFuncList2 f xs n =
   equating (take n) xs (zipWith (map . f) [0..] $ tail $ List.inits [0..])

binomials :: Int -> Bool
binomials = equalFuncList2 Comb.binomial Comb.binomials

catalanNumbers :: Int -> Bool
catalanNumbers = equalFuncList Comb.catalanNumber Comb.catalanNumbers

fibonacciNumbers :: Int -> Bool
fibonacciNumbers = equalFuncList Comb.fibonacciNumber Comb.fibonacciNumbers

derangementNumber :: Int -> Bool
derangementNumber = equalFuncList Comb.derangementNumber Comb.derangementNumbers

derangementNumbers :: Int -> Bool
derangementNumbers n =
   allEqual $ map (take n) $ (
      CombPriv.derangementNumbersPS0 :
      CombPriv.derangementNumbersPS1 :
      CombPriv.derangementNumbersInclExcl :
      [] :: [[Integer]])


bellSeries :: Int -> Bool
bellSeries =
   equalFuncList
      (\k -> round (Bell.bellSeries (fromInteger k) :: Double))
      (Bell.bellRec :: [Integer])


genPermutationWOFP :: QC.Gen (Int, String)
genPermutationWOFP = do
   xs <- take 6 . nub <$> QC.arbitrary
   k <- QC.choose (0, length xs)
   return (k,xs)

permutationWOFP :: Int -> String -> Bool
permutationWOFP k xs =
   PermWOFP.numbers !! length xs !! k == length (PermWOFP.enumerate k xs)

permutationWOFPFactorial :: Int -> Bool
permutationWOFPFactorial k =
   Comb.factorial (toInteger k) == PermWOFP.numbers !! k !! 0

permutationWOFPDerangement :: Int -> Bool
permutationWOFPDerangement k =
   Comb.derangementNumber (toInteger k) == PermWOFP.numbers !! k !! k


cardPairs1 :: Bool
cardPairs1 =
   case CardPairs.testCardsBorderDynamic of
      (x,y,z)  ->  x == y  &&  y == z

genCardCount :: QC.Gen (CardPairs.CardCount Int)
genCardCount =
   liftM3 CardPairs.CardCount
      (QC.choose (0,5)) (QC.choose (0,5)) (QC.choose (0,5))

cardPairs :: CardPairs.CardCount Int -> Bool
cardPairs cc =
   let x = CardPairs.possibilitiesCardsBorderNaive cc
       y = CardPairs.possibilitiesCardsBorderDynamic cc ! cc
       z = CardPairs.possibilitiesCardsBorder2Dynamic cc ! cc
   in  x == y  &&  y == z


genMastermindDistinct :: QC.Gen (Int, Int, Int, Int)
genMastermindDistinct = do
   n <- QC.choose (0,12)
   k <- QC.choose (0, min 5 n)
   b <- QC.choose (0,k)
   w <- QC.choose (0,k-b)
   return (n,k,b,w)

mastermindDistinct :: Int -> Int -> Int -> Int -> Bool
mastermindDistinct n k b w =
   let alphabet = take n ['a'..]
       code = take k alphabet
   in  Mastermind.numberDistinct n k b w ==
       (toInteger $ length $
        filter ((Mastermind.Eval b w ==) . Mastermind.evaluate code) $
        Comb.variate k alphabet)



testUnit :: Testable prop => String -> prop -> IO ()
testUnit label p = putStr (label++": ") >> quickCheck p

main :: IO ()
main =
   sequence_ $
      testUnit "permutation sums"
         (QC.forAll (take 6 <$> QC.arbitrary) permuteSum) :
      testUnit "permutations"
         (QC.forAll (take 6 <$> QC.arbitrary :: QC.Gen [Int]) permute) :
      testUnit "permuteAlt"
         (QC.forAll (take 6 <$> QC.arbitrary :: QC.Gen [Int]) permuteAlt) :
      testUnit "permuteRepM"
         (QC.forAll (genPermuteRep 10) permuteRepM) :
      testUnit "permuteRepNub"
         (QC.forAll (genPermuteRep  7) permuteRepNub) :
      testUnit "permuteRepNubBig"
         (QC.forAll (genPermuteRep 10) permuteRepNubBig) :
      testUnit "permuteRepMonotony"
         (QC.forAll (genPermuteRep 10) permuteRepMonotony) :
      testUnit "permuteRepChoose"
         (QC.forAll (QC.choose (0,10)) permuteRepChoose) :
      testUnit "choose"
         (QC.forAll genChoose (uncurry choose)) :
      testUnit "chooseLength"
         (QC.forAll (QC.choose (0,10)) chooseLength) :
      testUnit "chooseUnrank"
         (QC.forAll genChooseIndex $ uncurry3 chooseUnrank) :
      testUnit "chooseUnrankSequence"
         (QC.forAll (QC.choose (0,10)) chooseUnrankSequence) :
      testUnit "chooseRankUnrank"
         (QC.forAll genChooseIndex $ uncurry3 chooseRankUnrank) :
      testUnit "chooseUnrankRank" chooseUnrankRank :
      testUnit "variation without repetitions with list monad"
         (QC.forAll (QC.choose (-1,7)) $ \n ->
          QC.forAll genVariate $ variate n) :
      testUnit "variation with repetitions with list monad"
         (QC.forAll (QC.choose (-1,7)) $ \n ->
          QC.forAll genVariate $ variateRep n) :
      testUnit "variatePermute" (QC.forAll genVariate variatePermute) :
      testUnit "tuples" (QC.forAll genTuples (uncurry tuples)) :
      testUnit "permute expressed by variate"
         (variatePermuteClip 1000 :: String -> Bool) :
      testUnit "binomial vs. choose"
         (QC.forAll (QC.choose (0,12)) binomialChoose) :
      testUnit "multinomial vs. permutation with repetitions"
         (QC.forAll (genPermuteRep 10) multinomialPermuteRep) :
      testUnit "multinomial commutative"
         (QC.forAll (QC.listOf $ QC.choose (0,300)) multinomialCommutative) :
      testUnit "factorial vs. permute"
         (QC.forAll (take 8 <$> QC.arbitrary) factorial) :
      testUnit "binomial vs. tuples"
         (QC.forAll (take 16 <$> QC.arbitrary) binomial) :
      testUnit "binomial by factorial"
         (QC.forAll genBinomial $ uncurry binomialFactorial) :
      testUnit "factorial vs. factorials" (factorials 1000) :
      testUnit "binomial vs. binomials" (binomials 100) :
      testUnit "catalan numbers" (catalanNumbers 1000) :
      testUnit "fibonacci numbers" (fibonacciNumbers 10000) :
      testUnit "derangement number" (derangementNumber 1000) :
      testUnit "derangement numbers" (derangementNumbers 1000) :
      testUnit "set partition numbers"
         (QC.forAll (QC.choose (0,10000)) $ \n ->
          QC.forAll (take 7 <$> QC.arbitrary) $ setPartitionNumbers n) :
      testUnit "rectification numbers"
         (QC.forAll (QC.choose (0,7)) $ \n xs -> rectificationNumbers n xs) :
      testUnit "rectification montony"
         (QC.forAll (QC.choose (0,7)) $ \n xs ->
            rectificationsMonotony n (xs::[Int])) :
      testUnit "surjective mapping number" (surjectiveMappingNumber 20) :
      testUnit "surjective mapping numbers" (surjectiveMappingNumbers 20) :
      testUnit "bell series" (bellSeries 20) :
      testUnit "permutation without some fixpoints"
         (QC.forAll genPermutationWOFP $ uncurry permutationWOFP) :
      testUnit "permutation without some fixpoints vs. factorial"
         (QC.forAll (QC.choose (0,100)) permutationWOFPFactorial) :
      testUnit "permutation without some fixpoints vs. derangement"
         (QC.forAll (QC.choose (0,100)) permutationWOFPDerangement) :
      testUnit "partitions infinite linear factors"
         (QC.forAll (QC.choose (0,100)) Parts.propInfProdLinearFactors) :
      testUnit "partitions pentagonal power series"
         (Parts.propPentagonalPowerSeries 1000) :
      testUnit "partitions positive pentagonal numbers"
         (Parts.propPentagonalsDifP 10000) :
      testUnit "partitions negative pentagonal numbers"
         (Parts.propPentagonalsDifN 10000) :
      testUnit "partitions"
         (QC.forAll (QC.choose (1,10)) $ \k ->
          QC.forAll (QC.choose (0,50)) $ \n -> Parts.propPartitions k n) :
      testUnit "partitions count" (Parts.propNumPartitions 30) :
      testUnit "card pairs" cardPairs1 :
      testUnit "card pairs many" (QC.forAll genCardCount cardPairs) :
      testUnit "mastermind with distinct symbols"
         (QC.forAll genMastermindDistinct $ \(n,k,b,w) ->
            mastermindDistinct n k b w) :
      []
