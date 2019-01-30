{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE GADTs #-}
module Test.Banded (testsVar) where

import qualified Test.Generator as Gen
import qualified Test.Utility as Util
import Test.Banded.Utility
         (Square(Square), genSquare, genSquareCond,
          offDiagonals, offDiagonalNats)
import Test.Generator ((<.*|>), (<|*.>), (<.*.>), (<|*|>), (<|\|>))
import Test.Utility
         (approx, approxArray, approxMatrix,
          genOrder, genArray, Tagged, equalListWith)

import qualified Numeric.LAPACK.Matrix.Banded as Banded
import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix (ZeroInt, (<#>), (<#), (#>))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, absolute)

import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary.Literal as TypeNum
import qualified Type.Data.Num.Unary.Proof as Proof
import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num.Unary (unary, (:+:))

import qualified Data.Array.Comfort.Shape as Shape

import Foreign.Storable (Storable)

import Control.Applicative ((<$>))

import qualified Test.QuickCheck as QC


data Banded height width a =
   forall sub super.
   (Unary.Natural sub, Unary.Natural super) =>
   Banded (Banded.General sub super height width a)

instance
   (Show width, Show height, Show a,
    Shape.C width, Shape.C height, Storable a) =>
      Show (Banded height width a) where
   showsPrec p (Banded a) = showsPrec p a


genBanded ::
   (Class.Floating a) => Gen.Matrix a Int Int (Banded ZeroInt ZeroInt a)
genBanded =
      flip Gen.mapGenDim Gen.matrixDims $ \maxElem maxDim (height,width) -> do
   order <- genOrder
   kl <- QC.choose (0, toInteger maxDim)
   ku <- QC.choose (0, toInteger maxDim)
   Unary.reifyNatural kl $ \sub ->
      Unary.reifyNatural ku $ \super ->
      fmap Banded $ genArray maxElem $
         MatrixShape.bandedGeneral (unary sub, unary super) order height width

multiplyFullIdentity ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Banded ZeroInt ZeroInt a -> Bool
multiplyFullIdentity (Banded m) =
   let a = Banded.toFull m
   in approxArray a $
      Banded.multiplyFull m $ Square.toGeneral $ Square.identityFromWidth a


multiplyVectorDot ::
   (Class.Floating a, Eq a) =>
   (Vector ZeroInt a,
    Banded ZeroInt ZeroInt a,
    Vector ZeroInt a) ->
   Bool
multiplyVectorDot (x, Banded m, y) =
   Vector.dot x (m#>y) == Vector.dot (x<#m) y


multiplyFullAny ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Banded ZeroInt ZeroInt a,
    Matrix.General ZeroInt ZeroInt a) ->
   Bool
multiplyFullAny (Banded a, b) =
   approxArray
      (Banded.multiplyFull a b)
      (Matrix.multiply (Banded.toFull a) b)

multiplyFullColumns ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Banded ZeroInt ZeroInt a,
    Matrix.General ZeroInt ZeroInt a) ->
   Bool
multiplyFullColumns (Banded a, b) =
   equalListWith approxArray
      (Matrix.toColumns (Banded.multiplyFull a b))
      (map (Banded.multiplyVector a) (Matrix.toColumns b))


multiplyFullAssoc ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Banded ZeroInt ZeroInt a,
    Matrix.General ZeroInt ZeroInt a,
    Matrix.General ZeroInt ZeroInt a) ->
   Bool
multiplyFullAssoc (Banded a, b, c) =
   approxArray
      (Matrix.multiply (Banded.multiplyFull a b) c)
      (Banded.multiplyFull a (Matrix.multiply b c))


addOffDiagonals ::
   (Unary.Natural subA, Unary.Natural superA,
    Unary.Natural subB, Unary.Natural superB) =>
   Banded.General subA superA heightA widthA a ->
   Banded.General subB superB heightB widthB a ->
   (Proof.Nat (subA :+: subB), Proof.Nat (superA :+: superB))
addOffDiagonals a b =
   fst $ MatrixShape.addOffDiagonals (offDiagonals a) (offDiagonals b)

multiplyBanded ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Banded ZeroInt ZeroInt a,
    Banded ZeroInt ZeroInt a) ->
   Bool
multiplyBanded (Banded a, Banded b) =
   case addOffDiagonals a b of
      (Proof.Nat, Proof.Nat) ->
         approxArray
            (Banded.toFull (Banded.multiply a b))
            (Banded.multiplyFull a (Banded.toFull b))

multiplyBandedVectorAssoc ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Banded ZeroInt ZeroInt a,
    Banded ZeroInt ZeroInt a,
    Vector ZeroInt a) ->
   Bool
multiplyBandedVectorAssoc (Banded a, Banded b, x) =
   case addOffDiagonals a b of
      (Proof.Nat, Proof.Nat) ->
         approxArray (a #> b #> x) (Banded.multiply a b #> x)


multiplyBandedAssoc ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Banded ZeroInt ZeroInt a,
    Banded ZeroInt ZeroInt a,
    Banded ZeroInt ZeroInt a) ->
   Bool
multiplyBandedAssoc (Banded a, Banded b, Banded c) =
   let ab = Banded.multiply a b
       bc = Banded.multiply b c
       (subA,superA) = offDiagonalNats a
       (subB,superB) = offDiagonalNats b
       (subC,superC) = offDiagonalNats c
   in case (addOffDiagonals a b, addOffDiagonals b c) of
         ((Proof.Nat, Proof.Nat), (Proof.Nat, Proof.Nat)) ->
            case ((addOffDiagonals ab c, addOffDiagonals a bc),
                  (Proof.addAssoc subA subB subC,
                   Proof.addAssoc superA superB superC)) of
               (((Proof.Nat, Proof.Nat), (Proof.Nat, Proof.Nat)),
                (Proof.AddAssoc, Proof.AddAssoc)) ->
                  approxArray (Banded.multiply a bc) (Banded.multiply ab c)


data Upper size a =
   forall super. (Unary.Natural super) => Upper (Banded.Upper super size a)

instance
   (Show size, Show a, Shape.C size, Storable a) =>
      Show (Upper size a) where
   showsPrec p (Upper a) = showsPrec p a

genUpper :: (Class.Floating a) => Gen.Matrix a Int Int (Upper ZeroInt a)
genUpper = flip Gen.mapGenDim Gen.squareDim $ \maxElem maxDim size -> do
   order <- genOrder
   ku <- QC.choose (0, toInteger maxDim)
   Unary.reifyNatural ku $ \super ->
      fmap Upper $ genArray maxElem $
      MatrixShape.bandedSquare (unary TypeNum.u0, unary super) order size

multiplyUpperVector ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Upper ZeroInt a, Vector ZeroInt a) -> Bool
multiplyUpperVector (Upper m, x) =
   approxArray (m#>x) (Banded.toUpperTriangular m #> x)

multiplyLowerVector ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Upper ZeroInt a, Vector ZeroInt a) -> Bool
multiplyLowerVector (Upper up, x) =
   let lo = Banded.transpose up
   in approxArray (lo#>x) (Banded.toLowerTriangular lo #> x)


determinant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
determinant (Square a) =
   approx 0.5 (Banded.determinant a) (Square.determinant $ Banded.toFull a)


invertible ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
invertible (Square a) = absolute (Banded.determinant a) > 0.1

multiplySolve ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Square ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
multiplySolve (Square a, b) =
   approxMatrix 1e-2 (a <#> Banded.solve a b) b



checkForAll ::
   (Show a, QC.Testable test) =>
   Gen.T tag dim a -> (a -> test) -> Tagged tag QC.Property
checkForAll gen = Util.checkForAll (Gen.run gen 10 5)


testsVar ::
   (Show a, Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   ("multiplyFullIdentity",
      checkForAll genBanded multiplyFullIdentity) :
   ("multiplyFullAny",
      checkForAll ((,) <$> genBanded <|*|> Gen.matrix) multiplyFullAny) :
   ("multiplyVectorDot",
      checkForAll
         ((,,) <$> Gen.vector <.*|> genBanded <.*.> Gen.vector)
         multiplyVectorDot) :
   ("multiplyFullColumns",
      checkForAll ((,) <$> genBanded <|*|> Gen.matrix) multiplyFullColumns) :
   ("multiplyFullAssoc",
      checkForAll
         ((,,) <$> genBanded <|*|> Gen.matrix <|*|> Gen.matrix)
         multiplyFullAssoc) :
   ("multiplyBanded",
      checkForAll ((,) <$> genBanded <|*|> genBanded) multiplyBanded) :
   ("multiplyBandedVectorAssoc",
      checkForAll
         ((,,) <$> genBanded <|*|> genBanded <|*.> Gen.vector)
         multiplyBandedVectorAssoc) :
   ("multiplyBandedAssoc",
      checkForAll
         ((,,) <$> genBanded <|*|> genBanded <|*|> genBanded)
         multiplyBandedAssoc) :
   ("multiplyUpperVector",
      checkForAll ((,) <$> genUpper <|*.> Gen.vector) multiplyUpperVector) :
   ("multiplyLowerVector",
      checkForAll ((,) <$> genUpper <|*.> Gen.vector) multiplyLowerVector) :
   ("determinant",
      checkForAll genSquare determinant) :
   ("multiplySolve",
      checkForAll
         ((,) <$> genSquareCond invertible <|\|> Gen.matrix) multiplySolve) :
   []
