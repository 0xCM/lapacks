{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
module Test.Orthogonal (testsVar) where

import qualified Test.Generator as Gen
import qualified Test.Utility as Util
import Test.Generator ((<|*|>), (<|\|>))
import Test.Utility
         (approx, approxReal, approxArrayTol, approxMatrix, isIdentity, Tagged)

import qualified Numeric.LAPACK.Orthogonal.Householder as HH
import qualified Numeric.LAPACK.Orthogonal as Ortho
import qualified Numeric.LAPACK.Matrix.Hermitian as Herm
import qualified Numeric.LAPACK.Matrix.Triangular as Triangular
import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix.Square (Square)
import Numeric.LAPACK.Matrix (General, ZeroInt, (<#>))
import Numeric.LAPACK.Scalar (RealOf, absolute, selectReal)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)

import Control.Applicative (liftA2, (<$>))

import qualified Test.QuickCheck as QC


pseudoInverseProjection ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General ZeroInt ZeroInt a -> Bool
pseudoInverseProjection a =
   let ainv = snd $ Ortho.pseudoInverseRCond 1e-5 a
       tol = selectReal 1e-1 1e-5
   in approxArrayTol tol a (a <#> ainv <#> a) &&
      approxArrayTol tol ainv (ainv <#> a <#> ainv)

pseudoInverseHermitian ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General ZeroInt ZeroInt a -> Bool
pseudoInverseHermitian a =
   let ainv = snd $ Ortho.pseudoInverseRCond 1e-5 a
       tol = selectReal 1e-2 1e-5
       aainv = a <#> ainv
       ainva = ainv <#> a
   in approxMatrix tol aainv (Matrix.adjoint aainv) &&
      approxMatrix tol ainva (Matrix.adjoint ainva)

pseudoInverseFactored ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.Tall ZeroInt ZeroInt a,
    Matrix.Wide ZeroInt ZeroInt a) -> Bool
pseudoInverseFactored (a,b) =
   let pinv x = snd $ Ortho.pseudoInverseRCond 1e-5 x
   in approxMatrix (selectReal 1e-1 1e-5)
         (pinv (a <#> b)) (pinv b <#> pinv a)

pseudoInverseInverse ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
pseudoInverseInverse a =
   approxMatrix (selectReal 1e-1 1e-5)
      (Matrix.inverse a)
      (snd $ Ortho.pseudoInverseRCond 1e-5 a)


determinant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
determinant a =
   let detSquare = Square.determinant a
       detOrtho = Ortho.determinant a
   in approx
         (1e-3 * max 1 (max (absolute detSquare) (absolute detOrtho)))
         detSquare detOrtho

determinantAbsolute ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
determinantAbsolute a =
   let det = absolute $ Ortho.determinant a
       detAbs = Ortho.determinantAbsolute a
   in approxReal (1e-5 * max 1 (max det detAbs)) det detAbs

gramianDeterminant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General ZeroInt ZeroInt a -> Bool
gramianDeterminant a =
   let cov = Herm.covariance a
       Shape.ZeroBased n = Matrix.width a
       estimate = (Vector.sum (Herm.takeDiagonal cov) / fromIntegral n) ^ n
   in approxReal (1e-5 * max 1 estimate)
         (Herm.determinant cov)
         (Ortho.determinantAbsolute a ^ (2::Int))


multiplyDeterminantRight ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (General ZeroInt ZeroInt a, Square ZeroInt a) -> Bool
multiplyDeterminantRight (a,b) =
   let detA = Ortho.determinantAbsolute a
       detB = absolute $ Ortho.determinant b
   in approxReal
         (selectReal 1e-1 1e-5 * max 1 detA * max 1 detB)
         (Ortho.determinantAbsolute (a<#>b))
         (detA * detB)

multiplyDeterminantLeft ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Square ZeroInt a, General ZeroInt ZeroInt a) -> Bool
multiplyDeterminantLeft (a,b) =
   let detA = absolute $ Ortho.determinant a
       detB = Ortho.determinantAbsolute b
   in approxReal
         (selectReal 1e-1 1e-5 * max 1 detA * max 1 detB)
         (Ortho.determinantAbsolute (a<#>b))
         (detA * detB)


genFullRankTallRHS ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Gen.Matrix a Int Int
      (Matrix.Tall ZeroInt ZeroInt a,
       Matrix.General ZeroInt ZeroInt a)
genFullRankTallRHS = (,) <$> Gen.fullRankTall <|\|> Gen.matrix


normalEquationLeastSquares ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.Tall ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
normalEquationLeastSquares (a, b) =
   approxArrayTol
      (selectReal 10 1e-3)
      (Ortho.leastSquares a b)
      (Herm.solve (Herm.covariance $ Matrix.fromFull a) $
       Matrix.adjoint a <#> b)

specializedLeastSquares ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.Tall ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
specializedLeastSquares (a, b) =
   approxArrayTol
      (selectReal 1e-1 1e-5)
      (Ortho.leastSquares a b)
      (snd $ Ortho.leastSquaresMinimumNormRCond 1e-5 (Matrix.fromFull a) b)

householderLeastSquares ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.Tall ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
householderLeastSquares (a, b) =
   approxArrayTol
      (selectReal 1e-1 1e-5)
      (Ortho.leastSquares a b)
      (HH.leastSquares (HH.fromMatrix a) b)



genFullRankWideRHS ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Gen.Matrix a Int Int
      (Matrix.Wide ZeroInt ZeroInt a,
       Matrix.General ZeroInt ZeroInt a)
genFullRankWideRHS = (,) <$> Gen.fullRankWide <|\|> Gen.matrix


normalEquationMinimumNorm ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.Wide ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
normalEquationMinimumNorm (a, b) =
   approxArrayTol
      (selectReal 10 1e-3)
      (Ortho.minimumNorm a b)
      (Matrix.adjoint a <#>
       Herm.solve (Herm.covariance $ Matrix.fromFull $ Matrix.adjoint a) b)

specializedMinimumNorm ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.Wide ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
specializedMinimumNorm (a, b) =
   approxArrayTol
      (selectReal 1e-1 1e-5)
      (Ortho.minimumNorm a b)
      (snd $ Ortho.leastSquaresMinimumNormRCond 1e-5 (Matrix.fromFull a) b)

householderMinimumNorm ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.Wide ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
householderMinimumNorm (a, b) =
   approxArrayTol
      (selectReal 1e-1 1e-5)
      (Ortho.minimumNorm a b)
      (HH.minimumNorm (HH.fromMatrix $ Matrix.adjoint a) b)


complementDimension ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.Tall ZeroInt ZeroInt a -> Bool
complementDimension a =
   let b = Matrix.fromFull a Matrix.||| Matrix.fromFull (Ortho.complement a)
   in Shape.size (Matrix.height b) == Shape.size (Matrix.width b)

complementBiorthogonal ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.Tall ZeroInt ZeroInt a -> Bool
complementBiorthogonal a =
   all (approx 1e-3 0) $
   Array.toList $ Matrix.adjoint a <#> Ortho.complement a

complementOrthogonal ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.Tall ZeroInt ZeroInt a -> Bool
complementOrthogonal =
   isIdentity (selectReal 1e-3 1e-7) .
   Herm.toSquare . Herm.covariance . Matrix.fromFull . Ortho.complement


householderReconstruction ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.General ZeroInt ZeroInt a -> Bool
householderReconstruction a =
   approxArrayTol (selectReal 1e-3 1e-7)
      a (uncurry (<#>) (Ortho.householder a))

householderDeterminant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
householderDeterminant a =
   let detOrtho = Ortho.determinant a
       detHH = HH.determinant $ HH.fromMatrix a
   in approx 1e-5 detOrtho detHH


maybeConjugate ::
   (Shape.C sh, Class.Floating a) =>
   HH.Conjugation -> Array sh a -> Array sh a
maybeConjugate HH.NonConjugated = id
maybeConjugate HH.Conjugated = Vector.conjugate

maybeTranspose ::
   (Shape.C size, Class.Floating a, MatrixShape.TriDiag diag,
    MatrixShape.Content lo, MatrixShape.Content up) =>
   Herm.Transposition ->
   Triangular.Triangular up diag lo size a -> Square size a
maybeTranspose HH.NonTransposed = Triangular.toSquare
maybeTranspose HH.Transposed = Triangular.toSquare . Triangular.transpose

maybeAdjoint ::
   (Shape.C size, Class.Floating a) =>
   HH.Inversion -> Square size a -> Square size a
maybeAdjoint HH.NonInverted = id
maybeAdjoint HH.Inverted = Matrix.adjoint

householderSolveRR ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (HH.Transposition, HH.Conjugation) ->
   Matrix.Tall ZeroInt ZeroInt a -> Bool
householderSolveRR (trans,conj) a =
   let qr = HH.fromMatrix a
   in  isIdentity (selectReal 1e-3 1e-7) $
         HH.tallSolveR trans conj qr $
         maybeTranspose trans $ maybeConjugate conj $ HH.tallExtractR qr


householderMultiplyR ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   HH.Transposition ->
   (Matrix.Tall ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) ->
   Bool
householderMultiplyR trans (a,b) =
   let qr = HH.fromMatrix a
       r = HH.tallExtractR qr
   in approxArrayTol
         (selectReal 1e-3 1e-7)
         (HH.tallMultiplyR trans qr b)
         (case trans of
            HH.NonTransposed -> r <#> b
            HH.Transposed -> Triangular.transpose r <#> b)


householderQOrthogonal ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.General ZeroInt ZeroInt a -> Bool
householderQOrthogonal a =
   let q = HH.extractQ $ HH.fromMatrix a
   in isIdentity (selectReal 1e-3 1e-7) $ Matrix.adjoint q <#> q


householderMultiplyQ ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   HH.Inversion ->
   (Matrix.General ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) ->
   Bool
householderMultiplyQ inv (a,b) =
   let qr = HH.fromMatrix a
   in approxArrayTol
         (selectReal 1e-3 1e-7)
         (maybeAdjoint inv (HH.extractQ qr) <#> b)
         (HH.multiplyQ inv qr b)


householderTallQOrthogonal ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.Tall ZeroInt ZeroInt a -> Bool
householderTallQOrthogonal =
   isIdentity (selectReal 1e-3 1e-7) .
   Herm.toSquare . Herm.covariance . Matrix.fromFull .
   HH.tallExtractQ . HH.fromMatrix

householderTallMultiplyQ ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.Tall ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
householderTallMultiplyQ (a,b) =
   let qr = HH.fromMatrix a
   in approxArrayTol
         (selectReal 1e-3 1e-7)
         (HH.tallExtractQ qr <#> b)
         (HH.tallMultiplyQ qr b)

householderTallMultiplyQAdjoint ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.Tall ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
householderTallMultiplyQAdjoint (a,b) =
   let qr = HH.fromMatrix a
   in approxArrayTol
         (selectReal 1e-3 1e-7)
         (Matrix.adjoint (HH.tallExtractQ qr) <#> b)
         (HH.tallMultiplyQAdjoint qr b)



checkForAll ::
   (Show a, QC.Testable test) =>
   Gen.T tag dim a -> (a -> test) -> Tagged tag QC.Property
checkForAll gen = Util.checkForAll (Gen.run gen 3 5)


testsVar ::
   (Show a, Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   ("pseudoInverseProjection",
      checkForAll Gen.matrix pseudoInverseProjection) :
   ("pseudoInverseHermitian",
      checkForAll Gen.matrix pseudoInverseHermitian) :
   ("pseudoInverseFactored",
      checkForAll
         ((,) <$> Gen.fullRankTall <|*|> Gen.fullRankWide)
         pseudoInverseFactored) :
   ("pseudoInverseInverse",
      checkForAll Gen.invertible pseudoInverseInverse) :

   ("determinant",
      checkForAll Gen.square determinant) :
   ("determinantAbsolute",
      checkForAll Gen.square determinantAbsolute) :
   ("gramianDeterminant",
      checkForAll Gen.matrix gramianDeterminant) :
   ("multiplyDeterminantRight",
      checkForAll
         ((,) <$> Gen.matrix <|*|> Gen.square) multiplyDeterminantRight) :
   ("multiplyDeterminantLeft",
      checkForAll
         ((,) <$> (fst . Ortho.householder <$> Gen.square) <|*|> Gen.matrix)
         multiplyDeterminantLeft) :
   ("normalEquationLeastSquares",
      checkForAll genFullRankTallRHS normalEquationLeastSquares) :
   ("normalEquationMinimumNorm",
      checkForAll genFullRankWideRHS normalEquationMinimumNorm) :
   ("specializedLeastSquares",
      checkForAll genFullRankTallRHS specializedLeastSquares) :
   ("specializedMinimumNorm",
      checkForAll genFullRankWideRHS specializedMinimumNorm) :

   ("complementDimension",
      checkForAll Gen.tall complementDimension) :
   ("complementBiorthogonal",
      checkForAll Gen.tall complementBiorthogonal) :
   ("complementOrthogonal",
      checkForAll Gen.tall complementOrthogonal) :

   ("householderReconstruction",
      checkForAll Gen.matrix householderReconstruction) :
   ("householderDeterminant",
      checkForAll Gen.square householderDeterminant) :
   ("householderLeastSquares",
      checkForAll genFullRankTallRHS householderLeastSquares) :
   ("householderMinimumNorm",
      checkForAll genFullRankWideRHS householderMinimumNorm) :
   ("householderSolveRR",
      Gen.withExtra checkForAll
         (liftA2 (,) QC.arbitraryBoundedEnum QC.arbitraryBoundedEnum)
         Gen.fullRankTall householderSolveRR) :
   ("householderMultiplyR",
      Gen.withExtra checkForAll
         QC.arbitraryBoundedEnum ((,) <$> Gen.tall <|*|> Gen.matrix)
         householderMultiplyR) :
   ("householderQOrthogonal",
      checkForAll Gen.matrix householderQOrthogonal) :
   ("householderMultiplyQ",
      Gen.withExtra checkForAll
         QC.arbitraryBoundedEnum ((,) <$> Gen.matrix <|\|> Gen.matrix)
         householderMultiplyQ) :
   ("householderTallQOrthogonal",
      checkForAll Gen.tall householderTallQOrthogonal) :
   ("householderTallMultiplyQ",
      checkForAll ((,) <$> Gen.tall <|*|> Gen.matrix) householderTallMultiplyQ) :
   ("householderTallMultiplyQAdjoint",
      checkForAll
         ((,) <$> Gen.tall <|\|> Gen.matrix) householderTallMultiplyQAdjoint) :
   []
