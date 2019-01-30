module Combinatorics.Permutation.WithoutSomeFixpoints where

import Combinatorics (permute)

{- |
@enumerate n xs@ list all permutations of @xs@
where the first @n@ elements do not keep their position
(i.e. are no fixpoints).

This is a generalization of derangement.

Naive but comprehensible implementation.
-}
enumerate :: (Eq a) => Int -> [a] -> [[a]]
enumerate k xs = filter (and . zipWith (/=) xs . take k) $ permute xs

{- | <http://oeis.org/A047920> -}
numbers :: (Num a) => [[a]]
numbers =
   tail $ scanl (\row fac -> scanl (-) fac row) [] $
   scanl (*) 1 $ iterate (1+) 1
