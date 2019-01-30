{-# LANGUAGE ExistentialQuantification #-}
module Test.Banded.Utility where

import qualified Test.Generator as Gen
import Test.Utility (genOrder, genArray)

import qualified Numeric.LAPACK.Matrix.Banded as Banded
import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import Numeric.LAPACK.Matrix.Shape (UnaryProxy)
import Numeric.LAPACK.Matrix (ZeroInt)

import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary.Proof as Proof
import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num.Unary (unary)
import Type.Base.Proxy (Proxy(Proxy))

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape

import Foreign.Storable (Storable)

import Data.Tuple.HT (mapPair)

import qualified Test.QuickCheck as QC


-- cf. MatrixShape.Private
natFromProxy :: (Unary.Natural n) => UnaryProxy n -> Proof.Nat n
natFromProxy Proxy = Proof.Nat

offDiagonals ::
   Banded.Banded sub super vert horiz height width a ->
   (UnaryProxy sub, UnaryProxy super)
offDiagonals = MatrixShape.bandedOffDiagonals . Array.shape

offDiagonalNats ::
   (Unary.Natural sub, Unary.Natural super) =>
   Banded.Banded sub super vert horiz height width a ->
   (Proof.Nat sub, Proof.Nat super)
offDiagonalNats = mapPair (natFromProxy, natFromProxy) . offDiagonals


data Square size a =
   forall sub super.
   (Unary.Natural sub, Unary.Natural super) =>
   Square (Banded.Square sub super size a)

instance
   (Show size, Show a, Shape.C size, Storable a) =>
      Show (Square size a) where
   showsPrec p (Square a) = showsPrec p a

genSquare :: (Class.Floating a) => Gen.Matrix a Int Int (Square ZeroInt a)
genSquare = genSquareCond (const True)

genSquareCond ::
   (Class.Floating a) =>
   (Square ZeroInt a -> Bool) ->
   Gen.Matrix a Int Int (Square ZeroInt a)
genSquareCond cond =
      flip Gen.mapGenDim Gen.squareDim $ \maxElem maxDim size -> do
   order <- genOrder
   kl <- QC.choose (0, toInteger maxDim)
   ku <- QC.choose (0, toInteger maxDim)
   Unary.reifyNatural kl $ \sub ->
      Unary.reifyNatural ku $ \super ->
      (fmap Square $
         genArray maxElem $
            MatrixShape.bandedSquare (unary sub, unary super) order size)
      `QC.suchThat`
      cond
