module Numeric.LAPACK.Permutation (
   Permutation,
   Matrix.Inversion(..),
   fromPivots,
   toPivots,
   toMatrix,
   determinant,
   numberFromSign,
   transpose,
   multiply,
   apply,
   ) where

import Numeric.LAPACK.Permutation.Private
import qualified Numeric.LAPACK.Matrix.Private as Matrix
