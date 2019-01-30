{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE ExistentialQuantification #-}
{-# LANGUAGE GADTs #-}
module Test.BandedHermitian (testsVar) where

import qualified Test.Generator as Gen
import qualified Test.Utility as Util
import Test.Banded.Utility
         (Square(Square), genSquare, natFromProxy, offDiagonalNats)
import Test.Generator ((<.*|>), (<|*.>), (<.*.>), (<|*|>), (<|\|>))
import Test.Utility
         (approxReal, approxArray, approxRealArrayTol, approxMatrix,
          genOrder, genArray, Tagged, equalListWith)

import qualified Numeric.LAPACK.Matrix.BandedHermitianPositiveDefinite
                                                       as BandedHermitianPD
import qualified Numeric.LAPACK.Matrix.BandedHermitian as BandedHermitian
import qualified Numeric.LAPACK.Matrix.Banded as Banded
import qualified Numeric.LAPACK.Matrix.Hermitian as Hermitian
import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.ShapeStatic as ShapeStatic
import Numeric.LAPACK.Matrix (ZeroInt, zeroInt, (<#>), (<#), (#>))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, fromReal, absolute, selectReal)

import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary.Proof as Proof
import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num.Unary (unary)

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape

import Foreign.Storable (Storable)

import Control.Applicative (liftA2, (<$>))

import qualified Data.List.HT as ListHT
import Data.Traversable (for)
import Data.Tuple.HT (mapSnd)

import qualified Test.QuickCheck as QC


data BandedHermitian size a =
   forall offDiag.
   (Unary.Natural offDiag) =>
   BandedHermitian (BandedHermitian.BandedHermitian offDiag size a)

instance
   (Show size, Show a, Shape.C size, Storable a) =>
      Show (BandedHermitian size a) where
   showsPrec p (BandedHermitian a) = showsPrec p a


{-
Non-real elements on the diagonal.
-}
_genBandedHermitian ::
   (Class.Floating a) => Gen.Matrix a Int Int (BandedHermitian ZeroInt a)
_genBandedHermitian =
      flip Gen.mapGenDim Gen.squareDim $ \maxElem maxDim size -> do
   order <- genOrder
   k <- QC.choose (0, toInteger maxDim)
   Unary.reifyNatural k $ \numOff ->
      fmap BandedHermitian $ genArray maxElem $
         MatrixShape.bandedHermitian (unary numOff) order size

genBandedHermitian ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Gen.Matrix a Int Int (BandedHermitian ZeroInt a)
genBandedHermitian =
      flip Gen.mapGenDim Gen.squareDim $ \maxElem maxDim size -> do
   order <- genOrder
   k <- QC.choose (0, toInteger maxDim)
   Unary.reifyNatural k $ \numOff -> do
      let shape = MatrixShape.bandedHermitian (unary numOff) order size
      BandedHermitian . Array.fromList shape <$>
         (for (Shape.indices shape) $ \ix ->
            let real =
                  case ix of
                     MatrixShape.InsideBox r c -> r==c
                     MatrixShape.VertOutsideBox _ _ -> False
                     MatrixShape.HorizOutsideBox _ _ -> False
            in if real
                  then fromReal <$> Util.genReal maxElem
                  else Util.genElement maxElem)



convertToFull ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian ZeroInt a -> Bool
convertToFull (BandedHermitian a) =
   approxArray
      (Hermitian.toSquare $ BandedHermitian.toHermitian a)
      (Banded.toFull $ BandedHermitian.toBanded a)

takeDiagonal ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian ZeroInt a -> Bool
takeDiagonal (BandedHermitian a) =
   approxRealArrayTol 1e-5
      (Hermitian.takeDiagonal $ BandedHermitian.toHermitian a)
      (BandedHermitian.takeDiagonal a)

covariance ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
covariance (Square a) =
   let (sub,super) = offDiagonalNats a
   in case (Proof.addNat sub super, Proof.addComm sub super) of
         (Proof.Nat, Proof.AddComm) ->
            approxArray
               (BandedHermitian.toBanded $ BandedHermitian.covariance a)
               (Banded.adjoint a <#> a)



type StaticVector1 n = Vector (ShapeStatic.ZeroBased (Unary.Succ n))

data SumRank1 size a =
   forall offDiag.
   (Unary.Natural offDiag) =>
   SumRank1 size [(RealOf a, (Shape.Index size, StaticVector1 offDiag a))]

instance
   (Show size, Show (Shape.Index size), Show a, Show (RealOf a),
    Shape.C size, Storable a) =>
      Show (SumRank1 size a) where
   showsPrec p (SumRank1 sh a) = showsPrec p (sh,a)

genScaledVectors ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Gen.Vector a Int (SumRank1 ZeroInt a)
genScaledVectors =
   flip Gen.mapGen Gen.vectorDim $ \maxElem size@(Shape.ZeroBased n) -> do
      k <- QC.choose (0, n-1)
      Unary.reifyNatural (toInteger k) $ \numOff ->
         fmap (SumRank1 size) $
         if n==0
            then return []
            else
               QC.listOf $
                  liftA2 (,) (Util.genReal maxElem) $
                  liftA2 (,) (QC.choose (0,n-k-1))
                     (Util.genArray maxElem
                        (ShapeStatic.ZeroBased $ unary $ Unary.succ numOff))

sumRank1 ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order -> SumRank1 ZeroInt a -> Bool
sumRank1 order (SumRank1 sh xs) =
   approxArray
      (BandedHermitian.toHermitian $ BandedHermitian.sumRank1 order sh xs)
      (Hermitian.sumRank1 order sh $
       map (mapSnd (uncurry $ displace sh)) xs)

displace ::
   (Shape.C sh, Class.Floating a) =>
   ZeroInt -> Int -> Vector sh a -> Vector ZeroInt a
displace (Shape.ZeroBased n) k a =
   Array.mapShape (zeroInt . Shape.size) $
      Vector.constant (zeroInt k) 0
      `Vector.append`
      a
      `Vector.append`
      Vector.constant (zeroInt $ max 0 $ n - k - Shape.size (Array.shape a)) 0


multiplyIdentity ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian.Transposition -> Matrix.General ZeroInt ZeroInt a -> Bool
multiplyIdentity trans m =
   approxArray m
      (BandedHermitian.multiplyFull trans
         (BandedHermitian.identity (Matrix.height m)) m)

multiplyDiagonal ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian.Transposition ->
   (Vector ZeroInt ar, Matrix.General ZeroInt ZeroInt a) -> Bool
multiplyDiagonal trans (d,m) =
   approxArray
      (Matrix.scaleRowsReal d m)
      (BandedHermitian.multiplyFull trans (BandedHermitian.diagonal d) m)

multiplyFullIdentity ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian ZeroInt a -> Bool
multiplyFullIdentity (BandedHermitian m) =
   let a = Banded.toFull $ BandedHermitian.toBanded m
   in approxArray a $
      BandedHermitian.multiplyFull BandedHermitian.NonTransposed m $
      Square.identityFrom a


multiplyHermitianVector ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian.Transposition ->
   (BandedHermitian ZeroInt a, Vector ZeroInt a) ->
   Bool
multiplyHermitianVector trans (BandedHermitian m, x) =
   approxArray
      (BandedHermitian.multiplyVector trans m x)
      (Hermitian.multiplyVector trans (BandedHermitian.toHermitian m) x)

multiplyVectorDot ::
   (Class.Floating a, Eq a) =>
   (Vector ZeroInt a, BandedHermitian ZeroInt a, Vector ZeroInt a) -> Bool
multiplyVectorDot (x, BandedHermitian m, y) =
   Vector.dot x (m#>y) == Vector.dot (x<#m) y


multiplyFullAny ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian.Transposition ->
   (BandedHermitian ZeroInt a,
    Matrix.General ZeroInt ZeroInt a) ->
   Bool
multiplyFullAny trans (BandedHermitian a, b) =
   approxArray
      (BandedHermitian.multiplyFull trans a b)
      (Hermitian.multiplyFull trans (BandedHermitian.toHermitian a) b)

multiplyFullColumns ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian.Transposition ->
   (BandedHermitian ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
multiplyFullColumns trans (BandedHermitian a, b) =
   equalListWith approxArray
      (Matrix.toColumns (BandedHermitian.multiplyFull trans a b))
      (map (BandedHermitian.multiplyVector trans a) (Matrix.toColumns b))


multiplyFullAssoc ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian.Transposition ->
   (BandedHermitian ZeroInt a,
    Matrix.General ZeroInt ZeroInt a,
    Matrix.General ZeroInt ZeroInt a) ->
   Bool
multiplyFullAssoc trans (BandedHermitian a, b, c) =
   approxArray
      (Matrix.multiply (BandedHermitian.multiplyFull trans a b) c)
      (BandedHermitian.multiplyFull trans a (Matrix.multiply b c))



genBandedHPD ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Gen.Matrix a Int Int (BandedHermitian ZeroInt a)
genBandedHPD = flip Gen.mapGenDim Gen.squareDim $ \maxElem maxDim size -> do
   order <- genOrder
   kl <- QC.choose (0, toInteger maxDim)
   ku <- QC.choose (0, toInteger maxDim)
   Unary.reifyNatural kl $ \subU ->
      Unary.reifyNatural ku $ \superU ->
      let sub   = unary subU;   subP   = natFromProxy sub
          super = unary superU; superP = natFromProxy super
      in case (Proof.addNat subP superP, Proof.addComm subP superP) of
            (Proof.Nat, Proof.AddComm) ->
               fmap (BandedHermitian . BandedHermitian.covariance) $
                  (genArray maxElem $
                     MatrixShape.bandedSquare (sub, super) order size)
                  `QC.suchThat`
                  (\a -> absolute (Banded.determinant a) > 0.1)


determinant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian ZeroInt a -> Bool
determinant (BandedHermitian a) =
   let detB = BandedHermitianPD.determinant a
       detS = Hermitian.determinant $ BandedHermitian.toHermitian a
   in approxReal (selectReal 1 1e-3 * max 1 (abs detB + abs detS)) detB detS


multiplySolve ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (BandedHermitian ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
multiplySolve (BandedHermitian a, b) =
   approxMatrix (selectReal 10 1e-3) (a <#> BandedHermitianPD.solve a b) b

solveDecomposed ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (BandedHermitian ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
solveDecomposed (BandedHermitian a, b) =
   approxMatrix (selectReal 1e-3 1e-7)
      (BandedHermitianPD.solve a b)
      (BandedHermitianPD.solveDecomposed (BandedHermitianPD.decompose a) b)



eigenvaluesDeterminant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian ZeroInt a -> Bool
eigenvaluesDeterminant (BandedHermitian a) =
   let det = BandedHermitianPD.determinant a
       prod = Vector.product $ BandedHermitian.eigenvalues a
   in approxReal ((det+prod) * selectReal 0.5 1e-6) det prod

eigensystem ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian ZeroInt a -> Bool
eigensystem (BandedHermitian a) =
   let (q,d) = BandedHermitian.eigensystem a
   in  approxMatrix 1e-4
         (Banded.toFull $ BandedHermitian.toBanded a)
         (q <#> Matrix.scaleRowsReal d (Square.adjoint q))

eigenvaluesHermitian ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian ZeroInt a -> Bool
eigenvaluesHermitian (BandedHermitian a) =
   approxRealArrayTol (selectReal 1e-3 1e-5)
      (BandedHermitian.eigenvalues a)
      (Hermitian.eigenvalues $ BandedHermitian.toHermitian a)

eigensystemHermitian ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   BandedHermitian ZeroInt a -> QC.Property
eigensystemHermitian (BandedHermitian a) =
   let (q0,d0) = BandedHermitian.eigensystem a
       (q1,d1) = Hermitian.eigensystem $ BandedHermitian.toHermitian a
       unit = Matrix.adjoint q0 <#> q1
       tol = selectReal 1e-4 1e-7
   in not (or (ListHT.mapAdjacent (approxReal 0.1) (Array.toList d0)))
      QC.==>
      approxRealArrayTol tol d0 d1
      &&
      and
         (zipWith
            (\(r,c) x -> approxReal tol (absolute x) $ if r==c then 1 else 0)
            (Shape.indices $ Array.shape unit) (Array.toList unit))



checkForAll ::
   (Show a, QC.Testable test) =>
   Gen.T tag dim a -> (a -> test) -> Tagged tag QC.Property
checkForAll gen = Util.checkForAll (Gen.run gen 6 5)

checkForAllExtra ::
   (Show a, Show b, QC.Testable test) =>
   QC.Gen a -> Gen.T tag dim b ->
   (a -> b -> test) -> Tagged tag QC.Property
checkForAllExtra = Gen.withExtra checkForAll


testsVar ::
   (Show a, Show ar,
    Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   ("convertToFull",
      checkForAll genBandedHermitian convertToFull) :
   ("takeDiagonal",
      checkForAll genBandedHermitian takeDiagonal) :
   ("sumRank1",
      checkForAllExtra genOrder genScaledVectors sumRank1) :
   ("covariance",
      checkForAll genSquare covariance) :
   ("multiplyIdentity",
      checkForAllExtra QC.arbitraryBoundedEnum Gen.matrix multiplyIdentity) :
   ("multiplyDiagonal",
      checkForAllExtra QC.arbitraryBoundedEnum
         ((,) <$> Gen.vectorReal <.*|> Gen.matrix) multiplyDiagonal) :
   ("multiplyFullIdentity",
      checkForAll genBandedHermitian multiplyFullIdentity) :
   ("multiplyFullAny",
      checkForAllExtra QC.arbitraryBoundedEnum
         ((,) <$> genBandedHermitian <|*|> Gen.matrix) multiplyFullAny) :
   ("multiplyHermitianVector",
      checkForAllExtra QC.arbitraryBoundedEnum
         ((,) <$> genBandedHermitian <|*.> Gen.vector)
         multiplyHermitianVector) :
   ("multiplyVectorDot",
      checkForAll
         ((,,) <$> Gen.vector <.*|> genBandedHermitian <.*.> Gen.vector)
         multiplyVectorDot) :
   ("multiplyFullColumns",
      checkForAllExtra QC.arbitraryBoundedEnum
         ((,) <$> genBandedHermitian <|*|> Gen.matrix) multiplyFullColumns) :
   ("multiplyFullAssoc",
      checkForAllExtra QC.arbitraryBoundedEnum
         ((,,) <$> genBandedHermitian <|*|> Gen.matrix <|*|> Gen.matrix)
         multiplyFullAssoc) :

   ("determinant",
      checkForAll genBandedHPD determinant) :
   ("multiplySolve",
      checkForAll ((,) <$> genBandedHPD <|\|> Gen.matrix) multiplySolve) :
   ("solveDecomposed",
      checkForAll ((,) <$> genBandedHPD <|\|> Gen.matrix) solveDecomposed) :

   ("eigenvaluesDeterminant",
      checkForAll genBandedHPD eigenvaluesDeterminant) :
   ("eigensystem",
      checkForAll genBandedHermitian eigensystem) :
   ("eigenvaluesHermitian",
      checkForAll genBandedHermitian eigenvaluesHermitian) :
   ("eigensystemHermitian",
      checkForAll genBandedHermitian eigensystemHermitian) :
   []
