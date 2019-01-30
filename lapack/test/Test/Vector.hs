{-# LANGUAGE TypeFamilies #-}
module Test.Vector (testsVar) where

import qualified Test.Generator as Gen
import qualified Test.Utility as Util
import Test.Utility (Tagged(Tagged), TaggedGen)

import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.Scalar as Scalar
import Numeric.LAPACK.Matrix (ZeroInt, zeroInt)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf)

import qualified Numeric.Netlib.Class as Class

import Control.Applicative (liftA2, (<$>))

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable ((!))

import qualified Data.NonEmpty as NonEmpty
import Data.NonEmpty ((!:))

import qualified Test.QuickCheck as QC
import Test.ChasingBottoms.IsBottom (isBottom)


appendTakeDrop ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Int -> Vector ZeroInt a -> Bool
appendTakeDrop n x =
   Util.approxArray x $
   Array.mapShape (zeroInt . Shape.size)
      (Vector.append (Vector.take n x) (Vector.drop n x))

takeLeftRightAppend ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Vector ZeroInt a, Vector ZeroInt a) -> Bool
takeLeftRightAppend (x,y) =
   let xy = Vector.append x y
   in Util.approxArray x (Vector.takeLeft xy)
      &&
      Util.approxArray y (Vector.takeRight xy)


normInf ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Vector ZeroInt a -> Bool
normInf x =
   Vector.normInf x
   ==
   (NonEmpty.maximum $ 0 !: map Scalar.absolute (Array.toList x))

normInf1 ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Vector ZeroInt a -> Bool
normInf1 x =
   Vector.normInf1 x
   ==
   (NonEmpty.maximum $ 0 !: map Scalar.norm1 (Array.toList x))


genVector :: (Class.Floating a) => TaggedGen a (Vector ZeroInt a)
genVector = Tagged $ Util.genArray 10 . zeroInt =<< QC.choose (0,5)

normInfAppend ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar, RealOf ar ~ ar) =>
   (Vector ZeroInt a, Vector ZeroInt a) -> Bool
normInfAppend (x,y) =
   Vector.normInf (Vector.append x y)
   ==
   Vector.normInf (Vector.autoFromList [Vector.normInf x, Vector.normInf y])

normInf1Append ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar, RealOf ar ~ ar) =>
   (Vector ZeroInt a, Vector ZeroInt a) -> Bool
normInf1Append (x,y) =
   Vector.normInf1 (Vector.append x y)
   ==
   Vector.normInf1 (Vector.autoFromList [Vector.normInf1 x, Vector.normInf1 y])


argAbsMaximum ::
   (Eq a, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Vector ZeroInt a -> Bool
argAbsMaximum xs =
   let kx@(k,x) = Vector.argAbsMaximum xs
   in if Array.shape xs == zeroInt 0
         then isBottom kx
         else xs!k == x && Scalar.absolute x == Vector.normInf xs

argAbs1Maximum ::
   (Eq a, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Vector ZeroInt a -> Bool
argAbs1Maximum xs =
   let kx@(k,x) = Vector.argAbs1Maximum xs
   in if Array.shape xs == zeroInt 0
         then isBottom kx
         else xs!k == x && Scalar.norm1 x == Vector.normInf1 xs


checkForAll ::
   (Show a, QC.Testable test) =>
   Gen.T tag dim a -> (a -> test) -> Tagged tag QC.Property
checkForAll gen = Util.checkForAll (Gen.run gen 10 5)


testsVar ::
   (Show a,
    Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar, RealOf ar ~ ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   ("appendTakeDrop",
      Gen.withExtra checkForAll
         (QC.getNonNegative <$> QC.arbitrary) Gen.vector appendTakeDrop) :
   ("takeLeftRightAppend",
      Util.checkForAllPlain
         (liftA2 (liftA2 (,)) genVector genVector) takeLeftRightAppend) :
   ("normInf",
      checkForAll Gen.vector normInf) :
   ("normInf1",
      checkForAll Gen.vector normInf1) :
   ("normInfAppend",
      Util.checkForAllPlain
         (liftA2 (liftA2 (,)) genVector genVector) normInfAppend) :
   ("normInf1Append",
      Util.checkForAllPlain
         (liftA2 (liftA2 (,)) genVector genVector) normInf1Append) :
   ("argAbsMaximum",
      checkForAll Gen.vector argAbsMaximum) :
   ("argAbs1Maximum",
      checkForAll Gen.vector argAbs1Maximum) :
   []
