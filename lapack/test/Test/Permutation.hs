module Test.Permutation where

import qualified Numeric.LAPACK.Permutation as Perm
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Permutation (Inversion(Inverted, NonInverted))
import Numeric.LAPACK.Matrix (ZeroInt, zeroInt)
import Numeric.LAPACK.Vector (Vector)

import qualified Data.Array.Comfort.Storable as Array
import Data.Array.Comfort.Storable (Array)

import Foreign.C.Types (CInt)

import Control.Monad (forM)

import qualified Test.QuickCheck as QC


genPivots :: QC.Gen (Vector ZeroInt CInt)
genPivots = do
   nat <- QC.arbitrary
   let n = length nat
   let nc = fromIntegral n
   fmap (Vector.fromList (zeroInt n)) $
      forM (zip [1..] nat) $ \(i,()) -> QC.choose (i,nc)


permutationPivots :: Bool -> Array ZeroInt CInt -> Bool
permutationPivots dir xs =
   let inv = if dir then Inverted else NonInverted
   in Array.toList (Perm.toPivots inv (Perm.fromPivots inv (Array.shape xs) xs))
      ==
      Array.toList xs


tests :: [(String, QC.Property)]
tests =
   ("permutationPivots",
      QC.property $ QC.forAll genPivots . permutationPivots) :
   []
