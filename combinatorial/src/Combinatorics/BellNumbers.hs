module Combinatorics.BellNumbers where

import Combinatorics (binomials, )
import Combinatorics.Utility (scalarProduct, )
import qualified PowerSeries


{- List of Bell numbers computed with the recursive formula given in
   Wurzel 2004-06, page 136 -}
bellRec :: Num a => [a]
bellRec =
   1 : map (scalarProduct bellRec) binomials

bellSeries :: (Floating a, Enum a) => Int -> a
bellSeries n =
   scalarProduct
      (map (^n) [0..])
      (take 30 PowerSeries.derivativeCoefficients)
     / exp 1
