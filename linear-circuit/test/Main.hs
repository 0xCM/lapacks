module Main where

import qualified ResistorCube
import qualified Tree

import qualified Test.QuickCheck as QC


approx :: Double -> Double -> Bool
approx x y = abs (x-y) < 1e-8

test :: (QC.Testable prop) => String -> prop -> IO ()
test msg prop =
   putStr (msg ++ ": ") >> QC.quickCheck prop

main :: IO ()
main = do
   test "resistor cube" (approx ResistorCube.resistance (5/6))
   test "resistor tree"
      (\x -> approx (Tree.treeResistance x) (Tree.graphResistance x))
   test "orientation of resistors"
      (uncurry approx . Tree.flippedResistances)
