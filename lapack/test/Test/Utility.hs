{-# LANGUAGE TypeFamilies #-}
module Test.Utility where

import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix.Extent as Extent
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.Orthogonal as Ortho
import Numeric.LAPACK.Matrix.Square (Square)
import Numeric.LAPACK.Matrix.Shape (Order(RowMajor,ColumnMajor))
import Numeric.LAPACK.Matrix (ZeroInt)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, absolute)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)

import qualified Control.Monad.Trans.State as MS
import Control.Monad (replicateM)
import Control.Applicative (Applicative, liftA2, pure, (<*>), (<$>))

import qualified Data.List.HT as ListHT
import qualified Data.Complex as Complex
import Data.Complex (Complex((:+)))
import Data.Monoid (Monoid(mempty,mappend))
import Data.Semigroup (Semigroup((<>)))

import qualified Test.QuickCheck as QC
import Test.ChasingBottoms.IsBottom (isBottom)


equalListWith :: (a -> a -> Bool) -> [a] -> [a] -> Bool
equalListWith eq xs ys =
   and $ ListHT.takeWhileJust $
   zipWith
      (\mx my ->
         case (mx,my) of
            (Nothing,Nothing) -> Nothing
            (Just x, Just y) -> Just $ eq x y
            _ -> Just False)
      (map Just xs ++ repeat Nothing)
      (map Just ys ++ repeat Nothing)


approx ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) => ar -> a -> a -> Bool
approx tol x y = absolute (x-y) <= tol

approxReal :: (Class.Real a) => a -> a -> a -> Bool
approxReal tol x y = abs (x-y) <= tol


approxArrayTol ::
   (Shape.C shape, Eq shape, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   ar -> Array shape a -> Array shape a -> Bool
approxArrayTol tol x y =
   if Array.shape x == Array.shape y
     then and $ zipWith (approx tol) (Array.toList x) (Array.toList y)
     else error "approxArray: shapes mismatch"

approxArray ::
   (Shape.C shape, Eq shape, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Array shape a -> Array shape a -> Bool
approxArray = approxArrayTol 1e-5

approxRealArrayTol ::
   (Shape.C shape, Eq shape, Class.Real a) =>
   a -> Array shape a -> Array shape a -> Bool
approxRealArrayTol tol x y =
   if Array.shape x == Array.shape y
     then and $ zipWith (approxReal tol) (Array.toList x) (Array.toList y)
     else error "approxRealArray: shapes mismatch"

approxMatrix ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Eq width,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   ar ->
   Matrix.Full vert horiz height width a ->
   Matrix.Full vert horiz height width a -> Bool
approxMatrix tol x y =
   approxArrayTol tol
      (Matrix.toRowMajor $ Matrix.fromFull x)
      (Matrix.toRowMajor $ Matrix.fromFull y)


genReal :: (Class.Real a) => Integer -> QC.Gen a
genReal n = fromInteger <$> QC.choose (-n,n)

genComplex :: (Class.Real a) => Integer -> QC.Gen (Complex a)
genComplex n = liftA2 (Complex.:+) (genReal n) (genReal n)

genElement :: (Class.Floating a) => Integer -> QC.Gen a
genElement n =
   Class.switchFloating (genReal n) (genReal n) (genComplex n) (genComplex n)

genArray ::
   (Shape.C shape, Class.Floating a) =>
   Integer -> shape -> QC.Gen (Array shape a)
genArray n shape =
   Array.fromList shape <$> replicateM (Shape.size shape) (genElement n)


select :: [a] -> QC.Gen (a, [a])
select = QC.elements . ListHT.removeEach

genDistinct ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Integer -> Integer -> ZeroInt -> QC.Gen (Vector ZeroInt a)
genDistinct maxElemS maxElemD size@(Shape.ZeroBased n) = do
   let range k = map fromInteger [(-k)..k]
   xs <-
      MS.evalStateT (replicateM n $ MS.StateT select) $
      Class.switchFloating
         (range maxElemS)
         (range maxElemD)
         (liftA2 (:+) (range maxElemS) (range maxElemS))
         (liftA2 (:+) (range maxElemD) (range maxElemD))
   return $ Vector.fromList size xs


genOrder :: QC.Gen Order
genOrder = QC.elements [RowMajor, ColumnMajor]



invertible ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square sh a -> Bool
invertible a = absolute (Square.determinant a) > 0.1

fullRankTall ::
   (Shape.C height, Shape.C width,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.Tall height width a -> Bool
fullRankTall a = Ortho.determinantAbsolute a > 0.1


isIdentity ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   ar -> Square ZeroInt a -> Bool
isIdentity tol eye =
   approxArrayTol tol eye (Square.identityFrom eye)



newtype Tagged tag a = Tagged a
type TaggedGen tag a = Tagged tag (QC.Gen a)

instance Functor (Tagged tag) where
   fmap f (Tagged a) = Tagged (f a)

instance Applicative (Tagged tag) where
   pure = Tagged
   Tagged f <*> Tagged a = Tagged (f a)



checkForAllPlain ::
   (Show a, QC.Testable test) =>
   TaggedGen tag a -> (a -> test) -> Tagged tag QC.Property
checkForAllPlain (Tagged gen) test = Tagged $ QC.forAll gen test

checkForAll ::
   (Show a, QC.Testable test) =>
   TaggedGen tag (a, Match) -> (a -> test) -> Tagged tag QC.Property
checkForAll taggedGen test =
   checkForAllPlain taggedGen $ \(a,match) ->
      case match of
         Match -> QC.property $ test a
         Mismatch -> QC.property $ isBottom $ test a

{- |
In @DontForceMatch@ mode the test generators
may ignore generating matching dimensions.
If dimensions actually mismatch, a @Mismatch@ value is returned.
In this case the test driver asserts that
the test routine is aborted with an error.
However, a typical test type might be
\"generic implementation = specialized implementation\".
If the generic implementation correctly checks the sizes,
then the tester cannot detect a missing check in the specialized implementation.
So far the proposed way to avoid this problem
is to add a test that relies solely on the function to be tested.
If you have no better idea, compare an implementation with itself.
-}
data Match = Mismatch | Match
   deriving (Eq, Show)

instance Semigroup Match where
   (<>) = mappend

instance Monoid Match where
   mempty = Match
   mappend Match Match = Match
   mappend _ _ = Mismatch



prefix :: String -> [(String, test)] -> [(String, test)]
prefix msg =
   map (\(str,test) -> (msg ++ "." ++ str, test))
