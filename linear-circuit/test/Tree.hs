{- |
Arrange resistors according to a tree of parallel and serial compositions.
Compare resistance of trees with the general graph resistance computation.
-}
module Tree where

import qualified Math.LinearCircuit as LinearCircuit

import qualified Test.QuickCheck as QC

import qualified Data.Graph.Comfort as Graph
import Data.Graph.Comfort (Graph)

import qualified Control.Monad.Trans.Class as MT
import qualified Control.Monad.Trans.State as MS
import Control.Monad (liftM, liftM2, replicateM)

import qualified Data.Map as Map; import Data.Map (Map)

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty as NonEmpty
import qualified Data.Foldable as Fold
import qualified Data.List.Match as Match
import Data.Functor.Classes (Eq1, Ord1, Show1, eq1, compare1, showsPrec1)
import Data.Monoid (mappend)
import Data.Ord.HT (comparing)
import Data.Eq.HT (equating)


data T =
     Resistance Double
   | Serial (NonEmpty.T [] T)
   | Parallel (NonEmpty.T [] T)
   deriving (Show)

instance QC.Arbitrary T where
   arbitrary =
      let res = liftM Resistance $ QC.choose (0,1)
          go 0 = res
          go size =
            let subTree n =
                  let x = QC.resize (div size n) QC.arbitrary
                  in  liftM2 NonEmpty.cons x (replicateM (n-1) x)
            in  QC.frequency $
                  (3, res) :
                  (1, liftM Serial   (QC.choose (1,size) >>= subTree)) :
                  (1, liftM Parallel (QC.choose (1,size) >>= subTree)) :
                  []
      in  QC.sized go
   shrink tree =
      case tree of
         Resistance res ->
            let simpleRess = [0,1]
            in  if elem res simpleRess
                  then []
                  else map Resistance simpleRess
         Parallel xs -> NonEmpty.flatten xs ++ map Parallel (QC.shrink xs)
         Serial xs   -> NonEmpty.flatten xs ++ map Serial   (QC.shrink xs)


parallel2 :: Double -> Double -> Double
parallel2 0 0 = 0
parallel2 x y = x*y / (x+y)

treeResistance :: T -> Double
treeResistance x =
   case x of
      Resistance res -> res
      Serial xs -> Fold.foldl1 (+) $ fmap treeResistance xs
      Parallel xs -> Fold.foldl1 parallel2 $ fmap treeResistance xs


newtype EdgeId = EdgeId Int
   deriving (Eq, Ord, Show)

instance Enum EdgeId where
   fromEnum (EdgeId n) = n
   toEnum = EdgeId

newEdgeId :: (Monad m) => MS.StateT EdgeId m EdgeId
newEdgeId = do
   n <- MS.get
   MS.put $ succ n
   return n

data Edge a =
   Edge {
      edgeId :: EdgeId,
      edgeFrom, edgeTo :: a
   }
   deriving (Show)

instance Eq (Edge a) where (==) = equating edgeId
instance Ord (Edge a) where compare = comparing edgeId

instance Eq1 Edge where eq1 = (==)
instance Ord1 Edge where compare1 = compare
instance Show1 Edge where showsPrec1 = showsPrec

instance Fold.Foldable Edge where
   foldMap f (Edge _ x y) = mappend (f x) (f y)

instance Graph.Edge Edge where
   from (Edge _ n _) = n
   to (Edge _ _ n) = n

instance Graph.Reverse Edge where
   reverseEdge (Edge n from to) = Edge n to from


newtype Node = Node Int
   deriving (Eq, Ord, Show)

instance Enum Node where
   fromEnum (Node n) = n
   toEnum = Node

newNode :: (Monad m) => MS.StateT Node m Node
newNode = do
   n <- MS.get
   MS.put $ succ n
   return n

edgesFromTree ::
   T -> (Node, Node) ->
   MS.StateT EdgeId (MS.State Node) (Map (Edge Node) Double)
edgesFromTree tree (from, to) =
   case tree of
      Resistance res -> do
         e <- newEdgeId
         return $ Map.singleton (Edge e from to) res
      Serial xs -> do
         ns <- sequence $ Match.replicate (NonEmpty.tail xs) $ MT.lift newNode
         fmap Map.unions $ sequence $
            NonEmpty.flatten $
            NonEmptyC.zipWith edgesFromTree xs $
            NonEmpty.mapAdjacent (,) $
            NonEmpty.cons from $ NonEmpty.snoc ns to
      Parallel xs -> do
         fmap Map.unions $ mapM (flip edgesFromTree (from,to)) $
            NonEmpty.flatten xs

graphFromTree :: T -> (Graph Edge Node Double (), (Node, Node))
graphFromTree tree =
   let ((edgeMap, globalEnds), lastNode) =
         flip MS.runState (Node 0) $ flip MS.evalStateT (EdgeId 0) $ do
            ends <- MT.lift $ liftM2 (,) newNode newNode
            edges <- edgesFromTree tree ends
            return (edges, ends)
   in  (Graph.fromMap
          (Map.fromList $ map (flip (,) ()) [Node 0 .. pred lastNode])
          edgeMap,
        globalEnds)

graphResistance :: T -> Double
graphResistance =
   uncurry (uncurry . LinearCircuit.resistance) . graphFromTree



data
   FlippedGraph =
      FlippedGraph (Graph Edge Node (Double, Bool) ()) (Node, Node)
   deriving (Show)

instance QC.Arbitrary FlippedGraph where
   arbitrary = do
      (graph, ends) <- fmap graphFromTree QC.arbitrary
      flpGraph <-
         Graph.traverseEdge (\res -> liftM ((,) res) QC.arbitrary) graph
      return $ FlippedGraph flpGraph ends

flippedResistances :: FlippedGraph -> (Double, Double)
flippedResistances (FlippedGraph graph ends) =
   let flippedGraph =
          Graph.fromMap
             (Graph.nodeLabels graph)
             (Map.fromList $
              map
                 (\(e, (res, flp)) ->
                    (if flp then Graph.reverseEdge e else e, res)) $
              Map.toList $ Graph.edgeLabels graph)
   in  (uncurry (LinearCircuit.resistance (Graph.mapEdge fst graph)) ends,
        uncurry (LinearCircuit.resistance flippedGraph) ends)
