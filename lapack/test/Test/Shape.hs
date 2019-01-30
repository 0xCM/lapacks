{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ExistentialQuantification #-}
module Test.Shape where

import Test.Utility (genOrder, prefix)

import qualified Data.Array.Comfort.Shape.Test as ShapeTest
import qualified Data.Array.Comfort.Shape as Shape

import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent as Extent
import Numeric.LAPACK.Matrix (ZeroInt, zeroInt)

import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num.Unary (unary)

import Control.Applicative ((<$>))

import qualified Test.QuickCheck as QC


genGeneral :: QC.Gen (MatrixShape.General ZeroInt ZeroInt)
genGeneral = do
   order <- genOrder
   m <- QC.choose (0,10)
   n <- QC.choose (0,10)
   return $ MatrixShape.general order (zeroInt m) (zeroInt n)

genTall :: QC.Gen (MatrixShape.Tall ZeroInt ZeroInt)
genTall = do
   order <- genOrder
   m <- QC.choose (0,10)
   n <- QC.choose (0,m)
   return $ MatrixShape.tall order (zeroInt m) (zeroInt n)

genWide :: QC.Gen (MatrixShape.Wide ZeroInt ZeroInt)
genWide = do
   order <- genOrder
   m <- QC.choose (0,10)
   n <- QC.choose (m,10)
   return $ MatrixShape.wide order (zeroInt m) (zeroInt n)

genSquare :: QC.Gen (MatrixShape.Square ZeroInt)
genSquare = do
   order <- genOrder
   n <- QC.choose (0,10)
   return $ MatrixShape.square order (zeroInt n)


genHermitian :: QC.Gen (MatrixShape.Hermitian ZeroInt)
genHermitian = do
   order <- genOrder
   n <- QC.choose (0,10)
   return $ MatrixShape.hermitian order (zeroInt n)

genDiagonal :: QC.Gen (MatrixShape.Diagonal ZeroInt)
genDiagonal = do
   order <- genOrder
   n <- QC.choose (0,10)
   return $ MatrixShape.diagonal order (zeroInt n)

genLowerTriangular ::
   QC.Gen (MatrixShape.LowerTriangular MatrixShape.NonUnit ZeroInt)
genLowerTriangular = do
   order <- genOrder
   n <- QC.choose (0,10)
   return $ MatrixShape.lowerTriangular order (zeroInt n)

genUpperTriangular ::
   QC.Gen (MatrixShape.UpperTriangular MatrixShape.NonUnit ZeroInt)
genUpperTriangular = do
   order <- genOrder
   n <- QC.choose (0,10)
   return $ MatrixShape.upperTriangular order (zeroInt n)

genSymmetric :: QC.Gen (MatrixShape.Symmetric ZeroInt)
genSymmetric = do
   order <- genOrder
   n <- QC.choose (0,10)
   return $ MatrixShape.symmetric order (zeroInt n)


data Banded vert horiz height width =
   forall sub super.
   (Unary.Natural sub, Unary.Natural super) =>
   Banded (MatrixShape.Banded sub super vert horiz height width)

instance
   (Extent.C horiz, Extent.C vert,
    Show height, Show width, Shape.C height, Shape.C width) =>
      Show (Banded vert horiz height width) where
   showsPrec p (Banded sh) = showsPrec p sh

instance
   (Extent.C horiz, Extent.C vert, Shape.C height, Shape.C width) =>
      Shape.C (Banded vert horiz height width) where
   size (Banded sh) = Shape.size sh
   uncheckedSize (Banded sh) = Shape.uncheckedSize sh

instance
   (Extent.C horiz, Extent.C vert,
    Shape.Indexed height, Shape.Indexed width) =>
      Shape.Indexed (Banded vert horiz height width) where
   type Index (Banded vert horiz height width) =
            MatrixShape.BandedIndex (Shape.Index height) (Shape.Index width)
   indices (Banded sh) = Shape.indices sh
   offset (Banded sh) = Shape.offset sh
   uncheckedOffset (Banded sh) = Shape.uncheckedOffset sh
   inBounds (Banded sh) = Shape.inBounds sh

   sizeOffset (Banded sh) = Shape.sizeOffset sh
   uncheckedSizeOffset (Banded sh) = Shape.uncheckedSizeOffset sh

instance
   (Extent.C horiz, Extent.C vert,
    Shape.InvIndexed height, Shape.InvIndexed width) =>
      Shape.InvIndexed (Banded vert horiz height width) where

   indexFromOffset (Banded sh) = Shape.indexFromOffset sh
   uncheckedIndexFromOffset (Banded sh) = Shape.uncheckedIndexFromOffset sh


genBanded ::
   MatrixShape.Full vert horiz height width ->
   QC.Gen (Banded vert horiz height width)
genBanded sh = do
   kl <- QC.choose (0,10)
   ku <- QC.choose (0,10)
   Unary.reifyNatural kl $ \sub ->
      Unary.reifyNatural ku $ \super ->
      return $ Banded $ MatrixShape.bandedFromFull (unary sub, unary super) sh


data BandedHermitian size =
   forall offDiag.
   (Unary.Natural offDiag) =>
   BandedHermitian (MatrixShape.BandedHermitian offDiag size)

instance (Show size, Shape.C size) => Show (BandedHermitian size) where
   showsPrec p (BandedHermitian sh) = showsPrec p sh

instance (Shape.C size) => Shape.C (BandedHermitian size) where
   size (BandedHermitian sh) = Shape.size sh
   uncheckedSize (BandedHermitian sh) = Shape.uncheckedSize sh

instance (Shape.Indexed size) => Shape.Indexed (BandedHermitian size) where
   type Index (BandedHermitian size) =
            MatrixShape.BandedIndex (Shape.Index size) (Shape.Index size)
   indices (BandedHermitian sh) = Shape.indices sh
   offset (BandedHermitian sh) = Shape.offset sh
   uncheckedOffset (BandedHermitian sh) = Shape.uncheckedOffset sh
   inBounds (BandedHermitian sh) = Shape.inBounds sh

   sizeOffset (BandedHermitian sh) = Shape.sizeOffset sh
   uncheckedSizeOffset (BandedHermitian sh) = Shape.uncheckedSizeOffset sh

instance
   (Shape.InvIndexed size) => Shape.InvIndexed (BandedHermitian size) where

   indexFromOffset (BandedHermitian sh) =
      Shape.indexFromOffset sh
   uncheckedIndexFromOffset (BandedHermitian sh) =
      Shape.uncheckedIndexFromOffset sh


genBandedHermitian :: QC.Gen (BandedHermitian ZeroInt)
genBandedHermitian = do
   order <- genOrder
   n <- QC.choose (0,10)
   k <- QC.choose (0,10)
   Unary.reifyNatural k $ \numOff ->
      return $ BandedHermitian $
         MatrixShape.bandedHermitian (unary numOff) order (zeroInt n)


tests :: [(String, QC.Property)]
tests =
   prefix "General" (ShapeTest.tests genGeneral) ++
   prefix "Tall" (ShapeTest.tests genTall) ++
   prefix "Wide" (ShapeTest.tests genWide) ++
   prefix "Square" (ShapeTest.tests genSquare) ++

   prefix "Split.Reflector.General"
      (ShapeTest.tests $
       MatrixShape.splitFromFull MatrixShape.Reflector <$> genGeneral) ++
   prefix "Split.Reflector.Tall"
      (ShapeTest.tests $
       MatrixShape.splitFromFull MatrixShape.Reflector <$> genTall) ++
   prefix "Split.Reflector.Wide"
      (ShapeTest.tests $
       MatrixShape.splitFromFull MatrixShape.Reflector <$> genWide) ++
   prefix "Split.Reflector.Square"
      (ShapeTest.tests $
       MatrixShape.splitFromFull MatrixShape.Reflector <$> genSquare) ++
   prefix "Split.Triangle.General"
      (ShapeTest.tests $
       MatrixShape.splitFromFull MatrixShape.Triangle <$> genGeneral) ++
   prefix "Split.Triangle.Tall"
      (ShapeTest.tests $
       MatrixShape.splitFromFull MatrixShape.Triangle <$> genTall) ++
   prefix "Split.Triangle.Wide"
      (ShapeTest.tests $
       MatrixShape.splitFromFull MatrixShape.Triangle <$> genWide) ++
   prefix "Split.Triangle.Square"
      (ShapeTest.tests $
       MatrixShape.splitFromFull MatrixShape.Triangle <$> genSquare) ++

   prefix "Hermitian" (ShapeTest.tests genHermitian) ++
   prefix "Diagonal" (ShapeTest.tests genDiagonal) ++
   prefix "LowerTriangular" (ShapeTest.tests genLowerTriangular) ++
   prefix "UpperTriangular" (ShapeTest.tests genUpperTriangular) ++
   prefix "Symmetric" (ShapeTest.tests genSymmetric) ++

   prefix "Banded.General" (ShapeTest.tests $ genBanded =<< genGeneral) ++
   prefix "Banded.Tall" (ShapeTest.tests $ genBanded =<< genTall) ++
   prefix "Banded.Wide" (ShapeTest.tests $ genBanded =<< genWide) ++
   prefix "Banded.Square" (ShapeTest.tests $ genBanded =<< genSquare) ++
   prefix "BandedHermitian" (ShapeTest.tests genBandedHermitian) ++
   []
