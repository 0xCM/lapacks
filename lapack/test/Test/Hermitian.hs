{-# LANGUAGE TypeFamilies #-}
module Test.Hermitian (testsVar) where

import qualified Test.Generator as Gen
import qualified Test.Utility as Util
import Test.Generator ((<.*|>), (<|*.>), (<|*|>), (<|\|>))
import Test.Utility
         (approx, approxReal, approxArray, approxArrayTol, approxMatrix,
          Tagged, genOrder)

import qualified Numeric.LAPACK.Orthogonal.Householder as HH
import qualified Numeric.LAPACK.Matrix.HermitianPositiveDefinite as HermitianPD
import qualified Numeric.LAPACK.Matrix.Hermitian as Hermitian
import qualified Numeric.LAPACK.Matrix.Triangular as Triangular
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix.Hermitian (Hermitian)
import Numeric.LAPACK.Matrix.Square (Square)
import Numeric.LAPACK.Matrix.Shape (Order)
import Numeric.LAPACK.Matrix (General, ZeroInt, zeroInt, (<#), (<#>), (#>))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, fromReal, selectReal)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape

import Control.Applicative (liftA2, (<$>))

import qualified Data.NonEmpty.Class as NonEmptyC
import qualified Data.NonEmpty as NonEmpty

import qualified Test.QuickCheck as QC


covariance ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General ZeroInt ZeroInt a -> Bool
covariance x =
   approxArray
      (Matrix.fromFull $ Hermitian.toSquare $ Hermitian.covariance x)
      (Matrix.adjoint x <#> x)


outer ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Order -> Vector ZeroInt a -> Bool
outer order x =
   approxArray
      (Matrix.fromFull $ Hermitian.toSquare $ Hermitian.outer order x)
      (Matrix.outer order x x)


genScaledVectors ::
   (NonEmptyC.Gen f, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Gen.Vector a Int (ZeroInt, f (ar, Vector ZeroInt a))
genScaledVectors =
   flip Gen.mapGen Gen.vectorDim $ \maxElem size ->
      fmap ((,) size) $
      NonEmptyC.genOf $
         liftA2 (,) (Util.genReal maxElem) (Util.genArray maxElem size)

sumRank1 ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Order -> (ZeroInt, [(ar, Vector ZeroInt a)]) -> Bool
sumRank1 order (sh,xs) =
   approxArray
      (Matrix.fromFull $ Hermitian.toSquare $ Hermitian.sumRank1 order sh xs)
      (foldl Vector.add (Vector.constant (MatrixShape.general order sh sh) 0) $
       fmap (rank1 order) xs)

sumRank1NonEmpty ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Order -> NonEmpty.T [] (ar, Vector ZeroInt a) -> Bool
sumRank1NonEmpty order xs =
   approxArray
      (Matrix.fromFull $ Hermitian.toSquare $
       Hermitian.sumRank1NonEmpty order xs)
      (NonEmpty.foldl1 Vector.add $ fmap (rank1 order) xs)

rank1 ::
   (Eq size, Shape.C size, Class.Floating a) =>
   Order -> (RealOf a, Vector size a) -> Matrix.General size size a
rank1 order (r,x) = Vector.scaleReal r $ Matrix.outer order x x


genScaledVectorPairs ::
   (NonEmptyC.Gen f, Class.Floating a) =>
   Gen.Vector a Int (ZeroInt, f (a, (Vector ZeroInt a, Vector ZeroInt a)))
genScaledVectorPairs =
   flip Gen.mapGen Gen.vectorDim $ \maxElem size ->
      fmap ((,) size) $
      NonEmptyC.genOf $
         liftA2 (,) (Util.genElement maxElem) $
         liftA2 (,) (Util.genArray maxElem size) (Util.genArray maxElem size)

sumRank2 ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Order -> (ZeroInt, [(a, (Vector ZeroInt a, Vector ZeroInt a))]) -> Bool
sumRank2 order (sh,xys) =
   approxArray
      (Matrix.fromFull $ Hermitian.toSquare $ Hermitian.sumRank2 order sh xys)
      (foldl Vector.add (Vector.constant (MatrixShape.general order sh sh) 0) $
       fmap (rank2 order) xys)

sumRank2NonEmpty ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Order -> NonEmpty.T [] (a, (Vector ZeroInt a, Vector ZeroInt a)) -> Bool
sumRank2NonEmpty order xys =
   approxArray
      (Matrix.fromFull $ Hermitian.toSquare $
       Hermitian.sumRank2NonEmpty order xys)
      (NonEmpty.foldl1 Vector.add $ fmap (rank2 order) xys)

rank2 ::
   (Eq size, Shape.C size, Class.Floating a) =>
   Order -> (a, (Vector size a, Vector size a)) -> Matrix.General size size a
rank2 order (a,(x,y)) =
   let ax = Vector.scale a x
   in Vector.add
         (Matrix.outer order ax y)
         (Matrix.outer order y ax)


addAdjoint ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
addAdjoint x =
   approxArray
      (Hermitian.toSquare $ Hermitian.addAdjoint x)
      (Matrix.add (Matrix.adjoint x) x)


multiplySquare ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian ZeroInt a -> Bool
multiplySquare a =
   approxArray (Hermitian.toSquare $ Hermitian.square a) (a <#> a)

squareSquare ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian ZeroInt a -> Bool
squareSquare a =
   approxArray
      (Hermitian.toSquare $ Hermitian.square a)
      (Square.square $ Hermitian.toSquare a)

{-
multiplyPower ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Int, Hermitian ZeroInt a) -> Bool
multiplyPower (n,a) =
   let b = Hermitian.power (fromIntegral n) a
       c = nest n (Hermitian.multiply a) $ Hermitian.identityFrom a
   in approxArrayTol (1e-6 * (Vector.normInf1 b + Vector.normInf1 c)) b c
-}


multiplyVectorLeft ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Vector ZeroInt a, Hermitian ZeroInt a) -> Bool
multiplyVectorLeft (x,a) =
   approxArray (x <# Hermitian.toSquare a) (x <# a)

multiplyVectorRight ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Hermitian ZeroInt a, Vector ZeroInt a) -> Bool
multiplyVectorRight (a,x) =
   approxArray (Hermitian.toSquare a #> x) (a #> x)


multiplyLeft ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (General ZeroInt ZeroInt a, Hermitian ZeroInt a) -> Bool
multiplyLeft (a,b) =
   approxMatrix 1e-5 (a <#> Hermitian.toSquare b) (a <#> b)

multiplyRight ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Hermitian ZeroInt a, General ZeroInt ZeroInt a) -> Bool
multiplyRight (a,b) =
   approxArray (Hermitian.toSquare a <#> b) (a <#> b)


determinant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian ZeroInt a -> Bool
determinant a =
   approx
      (selectReal 1e-1 1e-5)
      (fromReal $ Hermitian.determinant a)
      (Square.determinant $ Hermitian.toSquare a)

choleskyQR ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.Tall ZeroInt ZeroInt a -> QC.Property
choleskyQR a =
   let qr = HH.fromMatrix a
       r = HH.tallExtractR qr
   in HH.determinantAbsolute qr > 0.1
      QC.==>
      approxArrayTol 1e-1
         (Matrix.scaleRows (Array.map signum $ Triangular.takeDiagonal r) $
          Triangular.toSquare r)
         (Triangular.toSquare $
          HermitianPD.decompose $ Hermitian.covariance $ Matrix.fromFull a)


invertible ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian sh a -> Bool
invertible a = abs (Hermitian.determinant a) > 0.1

inverse ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian ZeroInt a -> Bool
inverse a =
   approxArrayTol
      (selectReal 1 1e-5)
      (Hermitian.toSquare $ Hermitian.inverse a)
      (Square.inverse $ Hermitian.toSquare a)


solve ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Hermitian ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
solve (a, b) =
   approxMatrix (selectReal 1 1e-5)
      (Hermitian.solve a b)
      (Square.solve (Hermitian.toSquare a) b)



genPositiveDefinite ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Gen.Matrix a Int Int (Hermitian ZeroInt a)
genPositiveDefinite =
   flip Gen.mapGenDim Gen.squareDim $
         \maxElem maxDim width@(Shape.ZeroBased w) -> do
      height <- zeroInt <$> QC.choose (w,maxDim)
      order <- Util.genOrder
      Hermitian.covariance . Matrix.fromFull <$>
         Util.genArray maxElem (MatrixShape.tall order height width)
            `QC.suchThat` Util.fullRankTall

determinantPD ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian ZeroInt a -> Bool
determinantPD a =
   approxReal (selectReal 100 1e-4)
      (Hermitian.determinant a)
      (HermitianPD.determinant a)

inversePD ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian ZeroInt a -> Bool
inversePD a =
   approxArrayTol (selectReal 1000 1e-4)
      (Hermitian.inverse a)
      (HermitianPD.inverse a)

solvePD ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Hermitian ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
solvePD (a,b) =
   approxArrayTol (selectReal 1000 1e-4)
      (Hermitian.solve a b)
      (HermitianPD.solve a b)

solveDecomposedPD ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Hermitian ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
solveDecomposedPD (a,b) =
   approxArrayTol (selectReal 1e-1 1e-6)
      (HermitianPD.solve a b)
      (HermitianPD.solveDecomposed (HermitianPD.decompose a) b)



eigenvaluesDeterminant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian ZeroInt a -> Bool
eigenvaluesDeterminant a =
   approxReal
      (selectReal 1e-1 1e-5)
      (Hermitian.determinant a)
      (Vector.product $ Hermitian.eigenvalues a)

eigensystem ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian ZeroInt a -> Bool
eigensystem a =
   let (q,d) = Hermitian.eigensystem a
   in  approxMatrix 1e-4
         (Hermitian.toSquare a)
         (q <#> Matrix.scaleRowsReal d (Square.adjoint q))



checkForAll ::
   (Show a, QC.Testable test) =>
   Gen.T tag dim a -> (a -> test) -> Tagged tag QC.Property
checkForAll gen = Util.checkForAll (Gen.run gen 3 5)

checkForAllExtra ::
   (Show a, Show b, QC.Testable test) =>
   QC.Gen a -> Gen.T tag dim b ->
   (a -> b -> test) -> Tagged tag QC.Property
checkForAllExtra = Gen.withExtra checkForAll


testsVar ::
   (Show a, Show ar, Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   ("covariance",
      checkForAll Gen.matrix covariance) :
   ("outer",
      checkForAllExtra genOrder Gen.vector outer) :
   ("sumRank1",
      checkForAllExtra genOrder genScaledVectors sumRank1) :
   ("sumRank1NonEmpty",
      checkForAllExtra genOrder (snd <$> genScaledVectors) sumRank1NonEmpty) :
   ("sumRank2",
      checkForAllExtra genOrder genScaledVectorPairs sumRank2) :
   ("sumRank2NonEmpty",
      checkForAllExtra genOrder
         (snd <$> genScaledVectorPairs) sumRank2NonEmpty) :
   ("addAdjoint",
      checkForAll Gen.square addAdjoint) :
   ("multiplySquare",
      checkForAll Gen.hermitian multiplySquare) :
   ("squareSquare",
      checkForAll Gen.hermitian squareSquare) :

   ("multiplyVectorLeft",
      checkForAll ((,) <$> Gen.vector <.*|> Gen.hermitian) multiplyVectorLeft) :
   ("multiplyVectorRight",
      checkForAll ((,) <$> Gen.hermitian <|*.> Gen.vector) multiplyVectorRight) :
   ("multiplyLeft",
      checkForAll ((,) <$> Gen.matrix <|*|> Gen.hermitian) multiplyLeft) :
   ("multiplyRight",
      checkForAll ((,) <$> Gen.hermitian <|*|> Gen.matrix) multiplyRight) :

   ("determinant",
      checkForAll Gen.hermitian determinant) :
   ("choleskyQR",
      checkForAll Gen.tall choleskyQR) :

   ("inverse",
      checkForAll (Gen.hermitianCond invertible) inverse) :
   ("solve",
      checkForAll
         ((,) <$> Gen.hermitianCond invertible <|\|> Gen.matrix) solve) :

   ("determinantPD",
      checkForAll genPositiveDefinite determinantPD) :
   ("inversePD",
      checkForAll genPositiveDefinite inversePD) :
   ("solvePD",
      checkForAll ((,) <$> genPositiveDefinite <|\|> Gen.matrix) solvePD) :
   ("solveDecomposedPD",
      checkForAll
         ((,) <$> genPositiveDefinite <|\|> Gen.matrix) solveDecomposedPD) :

   ("eigenvaluesDeterminant",
      checkForAll Gen.hermitian eigenvaluesDeterminant) :
   ("eigensystem",
      checkForAll Gen.hermitian eigensystem) :
   []
