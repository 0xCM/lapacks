{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE Rank2Types #-}
{-# LANGUAGE EmptyDataDecls #-}
module Test.Generator where

import qualified Test.Logic as Logic
import qualified Test.Utility as Util
import Test.Logic (Dim, MatchMode(DontForceMatch,ForceMatch))
import Test.Utility (Match)

import qualified UniqueLogic.ST.TF.System.Simple as Sys

import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix.Triangular as Triangular
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix.Hermitian (Hermitian)
import Numeric.LAPACK.Matrix (ZeroInt, zeroInt)
import Numeric.LAPACK.Scalar (RealOf, fromReal, one)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)
import Data.Array.Comfort.Shape ((:+:))

import qualified Control.Monad.Trans.RWS as MRWS
import qualified Control.Applicative.HT as AppHT
import Control.Applicative (liftA2, (<*>), (<$>))

import Data.Traversable (for)
import Data.Tuple.HT (swap, mapFst3, mapThd3)

import qualified Test.QuickCheck as QC



{- |
@Cons generator@ with @generator maxElem@.
@generator@ constructs an array and maintains relations between the dimensions.
Dimensions will be choosen arbitrarily from the range @(0,maxDim)@.
Elements are choosen from the range @(-maxElem,maxElem)@.

I moved to the 's' tag to within the 'Cons' constructor
and furthermore defined 'TaggedVariables' to strip the 's' tag
from the Variables in 'dim'.
This way, we can easily define 'checkForAll' in the test modules.
Otherwise there would not be a way to quantify 'dim' while containing 's' tags.
That is, we would have to reset 'dim' to () before every call to 'checkForAll'.
-}
newtype T tag dim array =
   Cons {
      decons :: forall s.
         Integer ->
         Logic.M s (TaggedVariables s dim, Logic.System s, Logic.M s array)
   }

data Variable dim

type family TaggedVariables s tuple
type instance TaggedVariables s (Variable dim) = Logic.Variable s dim
type instance TaggedVariables s () = ()
type instance TaggedVariables s (a,b) =
                  (TaggedVariables s a, TaggedVariables s b)

instance Functor (T tag dim) where
   fmap f (Cons gen) = Cons $ \mapElem -> mapThd3 (fmap f) <$> gen mapElem

run ::
   T tag dim array -> Integer -> Int ->
   Util.TaggedGen tag (array, Match)
run gen maxElem maxDim =
   Util.Tagged $
      QC.elements [DontForceMatch, ForceMatch] >>=
      Logic.runSTInGen
         (do (_dim, sys, queries) <- decons gen maxElem
             Sys.solve sys
             queries)
         maxDim

withExtra ::
   (T tag dim (a,b) -> ((a,b) -> c) -> io) ->
   QC.Gen a -> T tag dim b -> (a -> b -> c) -> io
withExtra checkForAll genA genB test =
   checkForAll (mapGen (\_ b -> flip (,) b <$> genA) genB) (uncurry test)


mapGen ::
   (Integer -> a -> QC.Gen b) ->
   T tag dim a -> T tag dim b
mapGen f (Cons gen) =
   Cons $ \maxElem ->
      mapThd3 (Logic.liftGen . f maxElem =<<) <$> gen maxElem

mapGenDim ::
   (Integer -> Int -> a -> QC.Gen b) ->
   T tag dim a -> T tag dim b
mapGenDim f (Cons gen) =
   Cons $ \maxElem -> do
      (maxDim, _matchMode) <- Logic.M MRWS.ask
      mapThd3 (Logic.liftGen . f maxElem maxDim =<<) <$> gen maxElem


combine ::
   (forall s.
    TaggedVariables s dimF -> TaggedVariables s dimA ->
    (TaggedVariables s dimB, Logic.System s)) ->
   T tag dimF (a -> b) ->
   T tag dimA a ->
   T tag dimB b
combine combineDim (Cons genF) (Cons genA) =
   Cons $ \maxElem -> do
      (dimF,sysF,f) <- genF maxElem
      (dimA,sysA,a) <- genA maxElem
      let (dimB, constraint) = combineDim dimF dimA
      return (dimB, sysF >> sysA >> constraint, f <*> a)


type Scalar tag = T tag ()

scalar :: (Class.Floating a) => Scalar a a
scalar =
   Cons $ \ maxElem ->
      return ((), return (), Logic.liftGen $ Util.genElement maxElem)

(<.*.>) ::
   (Dim size, Eq size) =>
   Vector tag size (a -> b) ->
   Vector tag size a ->
   Scalar tag b
(<.*.>) = combine (\dimF dimA -> ((), Logic.ruleEqualDim dimF dimA))


queryZeroInt :: Logic.Variable s Int -> Logic.M s ZeroInt
queryZeroInt var = zeroInt <$> Logic.query var

type Vector tag size = T tag (Variable size)

vectorDim :: (Class.Floating a) => Vector a Int ZeroInt
vectorDim =
   Cons $ \ _maxElem -> do
      dim <- Sys.globalVariable
      return (dim, return (), queryZeroInt dim)

vector :: (Class.Floating a) => Vector a Int (Vector.Vector ZeroInt a)
vector = mapGen Util.genArray vectorDim

vectorReal ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Vector a Int (Vector.Vector ZeroInt ar)
vectorReal = mapGen Util.genArray vectorDim

(<.*|>) ::
   (Dim height, Eq height) =>
   Vector tag height (a -> b) ->
   Matrix tag height width a ->
   Vector tag width b
(<.*|>) =
   combine (\dim (height,width) -> (width, Logic.ruleEqualDim dim height))

(<|*.>) ::
   (Dim width, Eq width) =>
   Matrix tag height width (a -> b) ->
   Vector tag width a ->
   Vector tag height b
(<|*.>) =
   combine (\(height,width) dim -> (height, Logic.ruleEqualDim width dim))

(<.=.>) ::
   (Dim size, Eq size) =>
   Vector tag size (a -> b) ->
   Vector tag size a ->
   Vector tag size b
(<.=.>) =
   combine (\sizeF sizeA -> (sizeF, Logic.ruleEqualDim sizeF sizeA))


type Matrix tag height width = T tag (Variable height, Variable width)

matrixDims ::
   (Class.Floating a) => Matrix a Int Int (ZeroInt, ZeroInt)
matrixDims =
   Cons $ \ _maxElem -> do
      dims <- liftA2 (,) Sys.globalVariable Sys.globalVariable
      return (dims, return (), AppHT.mapPair (queryZeroInt,queryZeroInt) dims)

matrix ::
   (Class.Floating a) => Matrix a Int Int (Matrix.General ZeroInt ZeroInt a)
matrix =
   flip mapGen matrixDims $ \maxElem dims -> do
      order <- Util.genOrder
      Util.genArray maxElem $ uncurry (MatrixShape.general order) dims


squareDim :: (Class.Floating a) => Matrix a Int Int ZeroInt
squareDim =
   Cons $ \ _maxElem -> do
      dim <- Sys.globalVariable
      return ((dim,dim), return (), queryZeroInt dim)

squareShaped ::
   (Shape.C sh, Class.Floating a) =>
   (MatrixShape.Order -> ZeroInt -> sh) -> Matrix a Int Int (Array sh a)
squareShaped shape =
   flip mapGen squareDim $ \maxElem size -> do
      order <- Util.genOrder
      Util.genArray maxElem $ shape order size

square :: (Class.Floating a) => Matrix a Int Int (Square.Square ZeroInt a)
square = squareShaped MatrixShape.square

squareCond ::
   (Class.Floating a) =>
   (Square.Square ZeroInt a -> Bool) ->
   Matrix a Int Int (Square.Square ZeroInt a)
squareCond cond =
   flip mapGen squareDim $ \maxElem size -> do
      order <- Util.genOrder
      Util.genArray maxElem (MatrixShape.square order size)
         `QC.suchThat`
         cond

invertible ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix a Int Int (Square.Square ZeroInt a)
invertible = squareCond Util.invertible

diagonal ::
   (Class.Floating a) => Matrix a Int Int (Triangular.Diagonal ZeroInt a)
diagonal = squareShaped MatrixShape.diagonal

identity ::
   (MatrixShape.Content lo, MatrixShape.Content up, Class.Floating a) =>
   Matrix a Int Int (Triangular.Triangular lo MatrixShape.Unit up ZeroInt a)
identity =
   flip mapGen squareDim $ \ _maxElem size -> do
      order <- Util.genOrder
      return $ Triangular.identity order size

triangularCond ::
   (MatrixShape.Content up, MatrixShape.Content lo, MatrixShape.TriDiag diag,
    Class.Floating a) =>
   (Triangular.Triangular lo diag up ZeroInt a -> Bool) ->
   Matrix a Int Int (Triangular.Triangular lo diag up ZeroInt a)
triangularCond cond =
   flip mapGen squareDim $ \maxElem size -> do
      order <- Util.genOrder
      genTriangularArray maxElem
         (MatrixShape.Triangular
            MatrixShape.autoDiag MatrixShape.autoUplo order size)
         `QC.suchThat`
         cond

triangular ::
   (MatrixShape.Content up, MatrixShape.Content lo, MatrixShape.TriDiag diag,
    Class.Floating a) =>
   Matrix a Int Int (Triangular.Triangular lo diag up ZeroInt a)
triangular = triangularCond (const True)


newtype GenTriangularDiag lo up a diag =
   GenTriangularDiag {
      runGenTriangularDiag ::
         MatrixShape.Triangular lo diag up ZeroInt ->
         QC.Gen (Triangular.Triangular lo diag up ZeroInt a)
   }

genTriangularArray ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Class.Floating a) =>
   Integer ->
   MatrixShape.Triangular lo diag up ZeroInt ->
   QC.Gen (Triangular.Triangular lo diag up ZeroInt a)
genTriangularArray maxElem =
   runGenTriangularDiag $
   MatrixShape.switchTriDiag
      (GenTriangularDiag $ \shape ->
         Array.fromList shape <$>
            (for (Shape.indices shape) $ \(r,c) ->
               if r==c
                  then return one
                  else Util.genElement maxElem))
      (GenTriangularDiag $ Util.genArray maxElem)


tallDims :: (Class.Floating a) => Matrix a Int Int (ZeroInt, ZeroInt)
tallDims =
   Cons $ \ _maxElem -> do
      height <- Sys.globalVariable
      width  <- Sys.globalVariable
      return ((height,width), Logic.ruleLessOrEqual width height,
              liftA2 (,) (queryZeroInt height) (queryZeroInt width))

tall ::
   (Class.Floating a) =>
   Matrix a Int Int (Matrix.Tall ZeroInt ZeroInt a)
tall =
   flip mapGen tallDims $ \maxElem dims -> do
      order <- Util.genOrder
      Util.genArray maxElem $ uncurry (MatrixShape.tall order) dims

fullRankTall ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix a Int Int (Matrix.Tall ZeroInt ZeroInt a)
fullRankTall =
   flip mapGen tallDims $ \maxElem dims -> do
      order <- Util.genOrder
      Util.genArray maxElem (uncurry (MatrixShape.tall order) dims)
         `QC.suchThat` Util.fullRankTall


wide ::
   (Class.Floating a) =>
   Matrix a Int Int (Matrix.Wide ZeroInt ZeroInt a)
wide = Matrix.transpose <$> transpose tall

fullRankWide ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix a Int Int (Matrix.Wide ZeroInt ZeroInt a)
fullRankWide = Matrix.transpose <$> transpose fullRankTall


hermitian ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix a Int Int (Hermitian ZeroInt a)
hermitian = hermitianCond (const True)

hermitianCond ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Hermitian ZeroInt a -> Bool) ->
   Matrix a Int Int (Hermitian ZeroInt a)
hermitianCond cond =
   flip mapGen squareDim $ \maxElem size -> do
      order <- Util.genOrder
      let shape = MatrixShape.hermitian order size
      (Array.fromList shape <$>
         (for (Shape.indices shape) $ \(r,c) ->
            if r==c
               then fromReal <$> Util.genReal maxElem
               else Util.genElement maxElem))
         `QC.suchThat` cond


{-
There cannot be a pure/point function.
-}
(<|*|>) ::
   (Dim fuse, Eq fuse) =>
   Matrix tag height fuse (a -> b) ->
   Matrix tag fuse width a ->
   Matrix tag height width b
(<|*|>) =
   combine
      (\(height,fuseF) (fuseA,width) ->
         ((height,width), Logic.ruleEqualDim fuseF fuseA))

transpose ::
   Matrix tag height width a ->
   Matrix tag width height a
transpose (Cons gen) = Cons $ \maxElem -> mapFst3 swap <$> gen maxElem

(<|\|>) ::
   (Dim height, Eq height) =>
   Matrix tag height width (a -> b) ->
   Matrix tag height nrhs a ->
   Matrix tag width nrhs b
(<|\|>) a b = transpose a <|*|> b

(<***>) ::
   Vector tag height (a -> b) ->
   Vector tag width a ->
   Matrix tag height width b
(<***>) = combine (\height width -> ((height,width), return ()))


(<|=|>) ::
   (Dim height, Eq height) =>
   (Dim width, Eq width) =>
   Matrix tag height width (a -> b) ->
   Matrix tag height width a ->
   Matrix tag height width b
(<|=|>) =
   combine $ \(heightF,widthF) (heightA,widthA) ->
      ((heightF,widthF),
       Logic.ruleEqualDim heightF heightA >> Logic.ruleEqualDim widthF widthA)

(<===>) ::
   (Dim width, Eq width) =>
   Matrix tag heightA width (a -> b) ->
   Matrix tag heightB width a ->
   Matrix tag (heightA:+:heightB) width b
(<===>) (Cons genF) (Cons genA) =
   Cons $ \maxElem -> do
      heightB <- Sys.globalVariable
      ((heightF,widthF),sysF,f) <- genF maxElem
      ((heightA,widthA),sysA,a) <- genA maxElem
      return
         ((heightB,widthF),
          do sysF >> sysA
             Logic.ruleEqualDim widthF widthA
             Logic.ruleAppendDim heightF heightA heightB,
          f <*> a)

(<|||>) ::
   (Dim height, Eq height) =>
   Matrix tag height widthA (a -> b) ->
   Matrix tag height widthB a ->
   Matrix tag height (widthA:+:widthB) b
(<|||>) f a = transpose $ transpose f <===> transpose a


infixl 4 <.*.>, <.*|>, <|*.>, <|*|>, <|\|>, <***>, <.=.>, <|=|>, <===>, <|||>
