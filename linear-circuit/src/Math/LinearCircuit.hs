module Math.LinearCircuit (resistance) where

import qualified Data.Graph.Comfort as Graph
import Data.Graph.Comfort (Graph)

import qualified Numeric.Container as NC
import qualified Numeric.LinearAlgebra.HMatrix as HMatrix
import qualified Data.Packed.Matrix as Matrix
import qualified Data.Packed.Vector as Vector
import Numeric.LinearAlgebra.HMatrix (Field, (<\>))
import Data.Packed.Matrix (Matrix)
import Data.Packed.Vector (Vector)

import qualified Data.Map as Map
import qualified Data.List as List
import Data.Monoid (mconcat)

import Control.Functor.HT (outerProduct)
import Data.Bool.HT (if')


voltageMatrix ::
   (Graph.Edge edge, Ord node, Field a) =>
   Graph edge node a nodeLabel -> Matrix a
voltageMatrix gr =
   Matrix.fromLists $
   outerProduct
      (\e n ->
         if' (Graph.from e == n) 1 $
         if' (Graph.to   e == n) (-1) $
         0)
      (Graph.edges gr)
      (Graph.nodes gr)

{- |
It is almost currentMatrix = trans voltageMatrix,
except that a row is deleted in currentMatrix.
-}
currentMatrix ::
   (Graph.Edge edge, Ord node, Field a) =>
   Graph edge node a nodeLabel -> node -> node -> Matrix a
currentMatrix gr _src dst =
   Matrix.fromLists $
   outerProduct
      (\n e ->
         if' (Graph.from e == n) 1 $
         if' (Graph.to   e == n) (-1) $
         0)
      (List.delete dst $ Graph.nodes gr)
      (Graph.edges gr)

resistanceMatrix ::
   (Graph.Edge edge, Ord node, Field a) =>
   Graph edge node a nodeLabel -> Matrix a
resistanceMatrix gr =
   HMatrix.diag $ Vector.fromList $
   Map.elems $ Graph.edgeLabels gr

fullMatrix ::
   (Graph.Edge edge, Ord node, Field a) =>
   Graph edge node a nodeLabel -> node -> node -> Matrix a
fullMatrix gr src dst =
   let currents = currentMatrix gr src dst
       voltages = voltageMatrix gr
   in  Matrix.fromBlocks
          [[NC.konst 0 (1, Matrix.cols currents),
               Matrix.asRow $ Vector.fromList $
               map (\n -> if n==src then 1 else 0) $ Graph.nodes gr],
           [resistanceMatrix gr, voltages],
           [currents, NC.konst 0 (Matrix.rows currents, Matrix.cols voltages)]]

rhs ::
   (Graph.Edge edge, Ord node, Field a) =>
   Graph edge node a nodeLabel -> node -> node -> Vector a
rhs gr src dst =
   mconcat
      [NC.konst 0 1,
       NC.konst 0 (length (Graph.edges gr)),
       Vector.fromList $
       map (\n -> if n==src then 1 else 0) $
       List.delete dst $ Graph.nodes gr]


solution ::
   (Graph.Edge edge, Ord node, Field a) =>
   Graph edge node a nodeLabel -> node -> node -> Vector a
solution gr src dst =
   fullMatrix gr src dst <\> rhs gr src dst

resistance ::
   (Graph.Edge edge, Ord node, Field a) =>
   Graph edge node a nodeLabel -> node -> node -> a
resistance gr src dst =
   solution gr src dst
   `NC.atIndex`
   (length (Graph.edges gr) + length (takeWhile (dst/=) $ Graph.nodes gr))
