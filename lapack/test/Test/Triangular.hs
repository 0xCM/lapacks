{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE Rank2Types #-}
module Test.Triangular (testsVar) where

import qualified Test.Generator as Gen
import qualified Test.Utility as Util
import Test.Generator ((<.*|>), (<|*.>), (<|*|>), (<|\|>), (<|=|>))
import Test.Utility (approx, approxArray, approxArrayTol, approxMatrix, Tagged)

import qualified Numeric.LAPACK.Matrix.Triangular as Triangular
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix.Triangular (Triangular)
import Numeric.LAPACK.Matrix (General, ZeroInt, (<#), (<#>), (#>))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, selectReal, absolute)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable ((!))

import Control.Applicative ((<$>))

import Data.Traversable (for)
import Data.Tuple.HT (mapFst)

import qualified Test.QuickCheck as QC


forceOrder ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order ->
   (Triangular lo diag up ZeroInt a, Vector ZeroInt a) -> Bool
forceOrder order (a,x) =
   let ao = Triangular.forceOrder order a
   in MatrixShape.triangularOrder (Array.shape ao) == order
      &&
      approxArray (a #> x) (ao #> x)

addDistributive ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Eq lo, Eq up, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   ((Triangular lo diag up ZeroInt a,
     Triangular lo diag up ZeroInt a),
    Vector ZeroInt a) ->
   Bool
addDistributive ((a,b),x) =
   approxArray
      (Triangular.add
         (Triangular.strictNonUnitDiagonal a)
         (Triangular.strictNonUnitDiagonal b) #> x)
      (Vector.add (a#>x) (b#>x))

subDistributive ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Eq lo, Eq up, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   ((Triangular lo diag up ZeroInt a,
     Triangular lo diag up ZeroInt a),
    Vector ZeroInt a) ->
   Bool
subDistributive ((a,b),x) =
   approxArray
      (Triangular.sub
         (Triangular.strictNonUnitDiagonal a)
         (Triangular.strictNonUnitDiagonal b) #> x)
      (Vector.sub (a#>x) (b#>x))


multiplyIdentityVector ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, Vector ZeroInt a) -> Bool
multiplyIdentityVector (eye,a) =
   approxArray a (Triangular.multiplyVector eye a)

multiplyIdentityFull ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, General ZeroInt ZeroInt a) ->
   Bool
multiplyIdentityFull (eye,a) =
   approxArray a (Triangular.multiplyFull eye a)

multiplyIdentity ::
   (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
    Eq lo, Eq diag, Eq up,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, Triangular lo diag up ZeroInt a) ->
   Bool
multiplyIdentity (eye,a) =
   approxArray a (Triangular.multiply eye a)

multiplyVector ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, Vector ZeroInt a) -> Bool
multiplyVector (a,x) =
   approxArray
      (Triangular.toSquare a #> x)
      (Triangular.multiplyVector a x)

multiply ::
   (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, Triangular lo diag up ZeroInt a) ->
   Bool
multiply (a,b) =
   approxArray
      (Triangular.toSquare a <#> Triangular.toSquare b)
      (Triangular.toSquare $ Triangular.multiply a b)

multiplyFull ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, General ZeroInt ZeroInt a) ->
   Bool
multiplyFull (a,b) =
   approxArray
      (Triangular.toSquare a <#> b)
      (Triangular.multiplyFull a b)

multiplySquare ::
   (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Triangular lo diag up ZeroInt a -> Bool
multiplySquare a =
   approxArray
      (Triangular.toSquare $ Triangular.square a)
      (Triangular.multiplyFull a $ Triangular.toSquare a)

squareSquare ::
   (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Triangular lo diag up ZeroInt a -> Bool
squareSquare a =
   approxArray
      (Triangular.toSquare $ Triangular.square a)
      (Square.square $ Triangular.toSquare a)


multiplyVectorLeft ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Vector ZeroInt a, Triangular lo diag up ZeroInt a) -> Bool
multiplyVectorLeft (x,a) =
   approxArray (x <# Triangular.toSquare a) (x <# a)

multiplyVectorRight ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, Vector ZeroInt a) -> Bool
multiplyVectorRight (a,x) =
   approxArray (Triangular.toSquare a #> x) (a #> x)


multiplyLeft ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (General ZeroInt ZeroInt a, Triangular lo diag up ZeroInt a) -> Bool
multiplyLeft (a,b) =
   approxMatrix 1e-5 (a <#> Triangular.toSquare b) (a <#> b)

multiplyRight ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, General ZeroInt ZeroInt a) -> Bool
multiplyRight (a,b) =
   approxArray (Triangular.toSquare a <#> b) (a <#> b)



determinant ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Triangular lo diag up ZeroInt a -> Bool
determinant a =
   approx
      (selectReal 1e-1 1e-5)
      (Triangular.determinant a)
      (Square.determinant $ Triangular.toSquare a)


invertible ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Shape.C sh, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Triangular lo diag up sh a -> Bool
invertible a = absolute (Triangular.determinant a) > 0.1

genInvertible ::
   (MatrixShape.Content up, MatrixShape.Content lo, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   GenTriangular lo diag up a
genInvertible = Gen.triangularCond invertible

inverse ::
   (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Triangular lo diag up ZeroInt a -> Bool
inverse a =
   approxArrayTol
      (selectReal 1 1e-5)
      (Triangular.toSquare $ Triangular.inverse a)
      (Square.inverse $ Triangular.toSquare a)

inverseGeneric ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Triangular lo diag up ZeroInt a -> Bool
inverseGeneric a =
   approxArrayTol
      (selectReal 1 1e-5)
      (Triangular.toSquare $ Triangular.inverseGeneric a)
      (Square.inverse $ Triangular.toSquare a)


solve ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
solve (a, b) =
   approxMatrix (selectReal 1 1e-5)
      (Triangular.solve a b)
      (Square.solve (Triangular.toSquare a) b)

solveIdentity ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular lo diag up ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
solveIdentity (eye, a) =
   approxMatrix (selectReal 1e-3 1e-5)
      a (Triangular.solve eye a)



eigenvaluesDeterminant ::
   (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Triangular lo diag up ZeroInt a -> Bool
eigenvaluesDeterminant a =
   approx
      (selectReal 1e-1 1e-5)
      (Triangular.determinant a)
      (Vector.product $ Triangular.eigenvalues a)


genDiagonalizable ::
   (MatrixShape.Content lo, MatrixShape.Content up,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   GenTriangular lo MatrixShape.NonUnit up a
genDiagonalizable =
   flip Gen.mapGen Gen.squareDim $ \maxElem size -> do
      order <- Util.genOrder
      d <- Util.genDistinct 3 10 size
      let shape =
            MatrixShape.Triangular
               MatrixShape.NonUnit MatrixShape.autoUplo order size
      Array.fromList shape <$>
         (for (Shape.indices shape) $ \(r,c) ->
            if r==c
               then return (d!r)
               else Util.genElement maxElem)

eigensystem ::
   (MatrixShape.DiagUpLo lo up, Eq lo, Eq up,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order -> Triangular lo MatrixShape.NonUnit up ZeroInt a -> Bool
eigensystem order a =
   let (vr,d,vl) = Triangular.eigensystem a
       scal = Triangular.takeDiagonal $ Triangular.multiply vl vr
   in approxMatrix
         (selectReal 1e-3 1e-5)
         (Triangular.toSquare a)
         (Triangular.toSquare $
          vr
          `Triangular.multiply`
          Triangular.diagonal order (Vector.mul d $ Array.map recip scal)
          `Triangular.multiply`
          vl)


checkForAll ::
   (Show a, QC.Testable test) =>
   Gen.T tag dim a -> (a -> test) -> Tagged tag QC.Property
checkForAll gen = Util.checkForAll (Gen.run gen 3 5)

checkForAllExtra ::
   (Show a, Show b, QC.Testable test) =>
   QC.Gen a -> Gen.T tag dim b ->
   (a -> b -> test) -> Tagged tag QC.Property
checkForAllExtra = Gen.withExtra checkForAll


type GenTriangular lo diag up a =
      Gen.Matrix a Int Int (Triangular lo diag up ZeroInt a)


addSuperName :: String -> [(String, a)] -> [(String, a)]
addSuperName superName = map (mapFst ((superName++) . ("."++)))

checkAnyFlexDiag ::
   (MatrixShape.TriDiag diag) =>
   String ->
   (forall lo up.
    (MatrixShape.Content lo, MatrixShape.Content up,
     Eq lo, Eq up, Show lo, Show up) =>
    GenTriangular lo diag up a ->
    Tagged a QC.Property) ->
   (forall lo up.
    (MatrixShape.Content lo, MatrixShape.Content up,
     Eq lo, Eq up, Show lo, Show up) =>
    GenTriangular lo diag up a) ->
   [(String, Tagged a QC.Property)]
checkAnyFlexDiag name checker gen =
   (checkDiagUpLoFlexDiag name checker gen ++) $
   addSuperName name $
   ("Symmetric", checker (Triangular.asSymmetric <$> gen)) :
   []

checkAny ::
   String ->
   (forall lo up diag.
    (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
     Eq lo, Eq up, Show lo, Show up, Show diag) =>
    GenTriangular lo diag up a ->
    Tagged a QC.Property) ->
   (forall lo up diag.
    (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
     Eq lo, Eq up, Show lo, Show up, Show diag) =>
    GenTriangular lo diag up a) ->
   [(String, Tagged a QC.Property)]
checkAny name checker gen =
   checkAnyFlexDiag (name++".Unit") checker
      (Triangular.forceUnitDiagonal <$> gen) ++
   checkAnyFlexDiag (name++".NonUnit") checker
      (Triangular.forceNonUnitDiagonal <$> gen)


checkDiagUpLoFlexDiag ::
   (MatrixShape.TriDiag diag) =>
   String ->
   (forall lo up.
    (MatrixShape.DiagUpLo lo up, Eq lo, Eq up, Show lo, Show up) =>
    GenTriangular lo diag up a ->
    Tagged a QC.Property) ->
   (forall lo up.
    (MatrixShape.DiagUpLo lo up, Eq lo, Eq up, Show lo, Show up) =>
    GenTriangular lo diag up a) ->
   [(String, Tagged a QC.Property)]
checkDiagUpLoFlexDiag name checker gen =
   addSuperName name $
   ("Diagonal", checker (Triangular.asDiagonal <$> gen)) :
   ("Lower", checker (Triangular.asLower <$> gen)) :
   ("Upper", checker (Triangular.asUpper <$> gen)) :
   []

checkDiagUpLo ::
   String ->
   (forall lo up diag.
    (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
     Eq lo, Eq diag, Eq up, Show lo, Show diag, Show up) =>
    GenTriangular lo diag up a -> Tagged a QC.Property) ->
   (forall lo up diag.
    (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
     Eq lo, Eq diag, Eq up, Show lo, Show diag, Show up) =>
    GenTriangular lo diag up a) ->
   [(String, Tagged a QC.Property)]
checkDiagUpLo name checker gen =
   checkDiagUpLoFlexDiag (name++".Unit") checker
      (Triangular.forceUnitDiagonal <$> gen) ++
   checkDiagUpLoFlexDiag (name++".NonUnit") checker
      (Triangular.forceNonUnitDiagonal <$> gen)


testsVar ::
   (Show a, Show ar, Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   checkAny "forceOrder"
      (\gen ->
         checkForAllExtra Util.genOrder
            ((,) <$> gen <|*.> Gen.vector) forceOrder)
      Gen.triangular ++
   checkAny "addDistributive"
      (\gen ->
         checkForAll
            ((,) <$> ((,) <$> gen <|=|> gen) <|*.> Gen.vector)
            addDistributive)
      Gen.triangular ++
   checkAny "subDistributive"
      (\gen ->
         checkForAll
            ((,) <$> ((,) <$> gen <|=|> gen) <|*.> Gen.vector)
            subDistributive)
      Gen.triangular ++

   checkAny "multiplyIdentityVector"
      (\gen -> checkForAll ((,) <$> gen <|*.> Gen.vector) multiplyIdentityVector)
      (Triangular.relaxUnitDiagonal <$> Gen.identity) ++
   checkAny "multiplyIdentityFull"
      (\gen -> checkForAll ((,) <$> gen <|*|> Gen.matrix) multiplyIdentityFull)
      (Triangular.relaxUnitDiagonal <$> Gen.identity) ++
   checkDiagUpLo "multiplyIdentity"
      (\gen -> checkForAll ((,) <$> gen <|*|> Gen.triangular) multiplyIdentity)
      (Triangular.relaxUnitDiagonal <$> Gen.identity) ++
   checkAny "multiplyVector"
      (\gen -> checkForAll ((,) <$> gen <|*.> Gen.vector) multiplyVector)
      Gen.triangular ++
   checkAny "multiplyFull"
      (\gen -> checkForAll ((,) <$> gen <|*|> Gen.matrix) multiplyFull)
      Gen.triangular ++
   checkAny "multiplyVectorLeft"
      (\gen -> checkForAll ((,) <$> Gen.vector <.*|> gen) multiplyVectorLeft)
      Gen.triangular ++
   checkAny "multiplyVectorRight"
      (\gen -> checkForAll ((,) <$> gen <|*.> Gen.vector) multiplyVectorRight)
      Gen.triangular ++
   checkAny "multiplyLeft"
      (\gen -> checkForAll ((,) <$> Gen.matrix <|*|> gen) multiplyLeft)
      Gen.triangular ++
   checkAny "multiplyRight"
      (\gen -> checkForAll ((,) <$> gen <|*|> Gen.matrix) multiplyRight)
      Gen.triangular ++

   checkDiagUpLo "multiply"
      (\gen -> checkForAll ((,) <$> gen <|*|> gen) multiply)
      Gen.triangular ++
   checkDiagUpLo "multiplySquare"
      (\gen -> checkForAll gen multiplySquare)
      Gen.triangular ++
   checkDiagUpLo "squareSquare"
      (\gen -> checkForAll gen squareSquare)
      Gen.triangular ++

   checkAny "determinant"
      (\gen -> checkForAll gen determinant)
      Gen.triangular ++
   checkAny "solve"
      (\gen -> checkForAll ((,) <$> gen <|\|> Gen.matrix) solve)
      genInvertible ++
   checkAny "solveIdentity"
      (\gen -> checkForAll ((,) <$> gen <|\|> Gen.matrix) solveIdentity)
      (Triangular.relaxUnitDiagonal <$> Gen.identity) ++
   checkDiagUpLo "inverse"
      (\gen -> checkForAll gen inverse)
      genInvertible ++
   checkAny "inverseGeneric"
      (\gen -> checkForAll gen inverseGeneric)
      genInvertible ++

   checkDiagUpLo "eigenvaluesDeterminant"
      (\gen -> checkForAll gen eigenvaluesDeterminant)
      Gen.triangular ++
   checkDiagUpLoFlexDiag "eigensystem"
      (\gen -> checkForAllExtra Util.genOrder gen eigensystem)
      genDiagonalizable ++
   []
