module Numeric.LAPACK.Orthogonal.Householder (
   Householder,
   General,
   Tall,
   Wide,
   Square,
   mapExtent,
   fromMatrix,
   determinant,
   determinantAbsolute,
   leastSquares,
   minimumNorm,

   Matrix.Transposition(..),
   Matrix.Conjugation(..),
   Matrix.Inversion(..),
   extractQ,
   extractR,
   multiplyQ,

   tallExtractQ,
   tallExtractR,
   tallMultiplyQ,
   tallMultiplyQAdjoint,
   tallMultiplyR,
   tallSolveR,
   ) where

import qualified Numeric.LAPACK.Matrix.Private as Matrix
import Numeric.LAPACK.Orthogonal.Private
