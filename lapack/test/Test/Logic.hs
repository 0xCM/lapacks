{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE Rank2Types #-}
module Test.Logic where

import Test.Utility (Match(Match,Mismatch))

import qualified UniqueLogic.ST.TF.Rule as Rule
import qualified UniqueLogic.ST.TF.System.Simple as Sys
import qualified UniqueLogic.ST.TF.System as Sys (runApplyM)

import qualified Data.Ref as Ref
import Data.STRef (newSTRef, writeSTRef, readSTRef)

import Data.Array.Comfort.Shape ((:+:)((:+:)))

import qualified Control.Monad.Trans.Class as MT
import qualified Control.Monad.Trans.RWS as MRWS
import Control.Monad.ST (ST, runST)
import Control.Applicative (Applicative, liftA2, (<$>))

import qualified Test.QuickCheck as QC
import qualified Test.QuickCheck.GenT as GenT
import Test.QuickCheck.Gen (Gen(MkGen))
import Test.QuickCheck.GenT (GenT)

import System.Random (Random)


data MatchMode = DontForceMatch | ForceMatch
   deriving (Eq, Show)

newtype M s a = M {runM :: MRWS.RWST (Int, MatchMode) Match () (GenT (ST s)) a}
   deriving (Functor, Applicative, Monad)

instance Ref.C (M s) where
   new =
      Ref.newCons
         (liftST . newSTRef)
         ((liftST .) . writeSTRef)
         (liftST . readSTRef)

liftST :: ST s a -> M s a
liftST = M . MT.lift . MT.lift

liftGen :: QC.Gen a -> M s a
liftGen = M . MT.lift . GenT.liftGen


type Variable s a = Sys.Variable (M s) a
type System s = Sys.T (M s) ()


example :: Int -> MatchMode -> QC.Gen ([Int], Match)
example =
   runSTInGen (do
      a <- Sys.globalVariable
      b <- Sys.globalVariable
      c <- Sys.globalVariable
      d <- Sys.globalVariable
      e <- Sys.globalVariable
      f <- Sys.globalVariable
      Sys.solve $ do
         ruleLessOrEqual a b
         ruleEqualDim d f
         ruleEqualDim c d
      mapM query [a,b,c,d,e,f])


choose :: (Random a) => (Int -> (a,a)) -> M s a
choose f = liftGen . QC.choose . f . fst =<< M MRWS.ask

ruleLessOrEqual :: Variable s Int -> Variable s Int -> System s
ruleLessOrEqual va vb = do
   assignmentM (\a -> choose (\ maxk -> (a,maxk))) va vb
   assignmentM (\b -> choose (\_maxk -> (0,b))) vb va


class Dim dim where chooseDim :: M s dim
instance Dim Int where chooseDim = choose ((,) 0)
instance (Dim dimA, Dim dimB) => Dim (dimA:+:dimB) where
   chooseDim = liftA2 (:+:) chooseDim chooseDim

ruleEqualDim ::
   (Dim dim, Eq dim) => Variable s dim -> Variable s dim -> System s
ruleEqualDim va vb = do
   let equalM =
         assignmentM $ \x -> do
            matchMode <- M $ MRWS.asks snd
            case matchMode of
               ForceMatch -> return x
               DontForceMatch -> do
                  y <- chooseDim
                  M $ MRWS.tell $ if x==y then Match else Mismatch
                  return y
   equalM va vb
   equalM vb va

assignmentM ::
   (Ref.C s) => (a -> s b) -> Sys.Variable s a -> Sys.Variable s b -> Sys.T s ()
assignmentM f vx vy = Sys.runApplyM (f <$> Sys.arg vx) vy


ruleAppendDim ::
   Variable s dimA -> Variable s dimB -> Variable s (dimA:+:dimB) -> System s
ruleAppendDim va vb vab = do
   Sys.assignment3 (:+:) va vb vab
   Sys.assignment2 (\(a:+:_) -> a) vab va
   Sys.assignment2 (\(_:+:b) -> b) vab vb


runSTInGen :: (forall s. M s b) -> Int -> MatchMode -> QC.Gen (b, Match)
runSTInGen m =
   \maxDim matchMode -> MkGen $ \r n ->
      runST (GenT.unGenT (MRWS.evalRWST (runM m) (maxDim,matchMode) ()) r n)

query :: (Dim dim) => Variable s dim -> M s dim
query v = do
   mk <- Sys.query v
   case mk of
      Just k -> return k
      Nothing -> do
         k <- chooseDim
         Sys.solve $ Rule.equ v =<< Sys.constant k
         return k
