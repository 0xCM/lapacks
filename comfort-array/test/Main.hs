{-# LANGUAGE TypeFamilies #-}
module Main where

import qualified Test.Shape as TestShape
import Test.Utility (prefix)

import qualified Test.QuickCheck as QC


main :: IO ()
main =
   mapM_ (\(name,prop) -> putStr (name ++ ": ") >> QC.quickCheck prop) $
   prefix "Shape" TestShape.tests ++
   []
