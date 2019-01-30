{- |
Consider a cube of resistors of equal resistance.
What is the overall resistance from one corner to the opposite one?
-}
module ResistorCube where

import qualified Math.LinearCircuit as LinearCircuit

import qualified Data.Graph.Comfort as Graph
import Data.Graph.Comfort (Graph)

import Control.Applicative (liftA2, liftA3)


data Coord = C0 | C1 deriving (Eq, Ord, Show, Enum, Bounded)
data Corner = Corner Coord Coord Coord deriving (Eq, Ord, Show)


allCoords :: [Coord]
allCoords = [minBound .. maxBound]

dimEdges ::
   (Coord -> Coord -> Coord -> Corner) ->
   [(Graph.UndirEdge Corner, Double)]
dimEdges corner =
   liftA2
      (\a b -> (Graph.undirEdge (corner C0 a b) (corner C1 a b), 1))
      allCoords allCoords

graph :: Graph Graph.UndirEdge Corner Double ()
graph =
   Graph.fromList
      (map (flip (,) ()) $ liftA3 Corner allCoords allCoords allCoords)
      (dimEdges (\x y z -> Corner x y z) ++
       dimEdges (\y z x -> Corner x y z) ++
       dimEdges (\z x y -> Corner x y z))

resistance :: Double
resistance =
   LinearCircuit.resistance graph
      (Corner C0 C0 C0) (Corner C1 C1 C1)
