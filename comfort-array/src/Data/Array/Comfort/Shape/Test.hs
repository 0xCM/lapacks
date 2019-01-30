{-# LANGUAGE TypeFamilies #-}
module Data.Array.Comfort.Shape.Test (tests) where

import qualified Data.Array.Comfort.Shape as Shape
import Data.Tuple.HT (mapSnd)

import qualified Test.QuickCheck as QC


uncheckedSize :: (Shape.C sh) => sh -> Bool
uncheckedSize sh  =  Shape.size sh == Shape.uncheckedSize sh

inBounds :: (Shape.Indexed sh) => sh -> Bool
inBounds sh  =  all (Shape.inBounds sh) $ Shape.indices sh


forAllIndices ::
   (Shape.Indexed sh, Shape.Index sh ~ ix, Show ix, QC.Testable prop) =>
   sh -> (ix -> prop) -> QC.Property
forAllIndices sh f =
   let ixs = Shape.indices sh
   in not (null ixs)  QC.==>  QC.forAll (QC.elements ixs) f

sizeOffset ::
   (Shape.Indexed sh, Shape.Index sh ~ ix, Show ix) => sh -> QC.Property
sizeOffset sh =
   forAllIndices sh $ \ix ->
      mapSnd ($ix) (Shape.sizeOffset sh)
      ==
      (Shape.size sh, Shape.offset sh ix)

uncheckedSizeOffset ::
   (Shape.Indexed sh, Shape.Index sh ~ ix, Show ix) => sh -> QC.Property
uncheckedSizeOffset sh =
   forAllIndices sh $ \ix ->
      mapSnd ($ix) (Shape.uncheckedSizeOffset sh) ==
         (Shape.uncheckedSize sh, Shape.uncheckedOffset sh ix)

uncheckedOffset ::
   (Shape.Indexed sh, Shape.Index sh ~ ix, Show ix) => sh -> QC.Property
uncheckedOffset sh =
   forAllIndices sh $ \ix ->
      Shape.offset sh ix == Shape.uncheckedOffset sh ix

lengthIndices :: (Shape.Indexed sh) => sh -> Bool
lengthIndices sh  =  length (Shape.indices sh) == Shape.size sh

indexOffsets :: (Shape.Indexed sh) => sh -> Bool
indexOffsets sh =
   map (Shape.offset sh) (Shape.indices sh) == take (Shape.size sh) [0..]

invIndices :: (Shape.InvIndexed sh, Shape.Index sh ~ ix, Eq ix) => sh -> Bool
invIndices sh =
   Shape.indices sh ==
   map (Shape.indexFromOffset sh) (take (Shape.size sh) [0..])

uncheckedInvIndices ::
   (Shape.InvIndexed sh, Shape.Index sh ~ ix, Eq ix) => sh -> Bool
uncheckedInvIndices sh =
   Shape.indices sh ==
   map (Shape.uncheckedIndexFromOffset sh) (take (Shape.size sh) [0..])


tests ::
   (Shape.InvIndexed sh, Show sh, Shape.Index sh ~ ix, Eq ix, Show ix) =>
   QC.Gen sh -> [(String, QC.Property)]
tests gen =
   ("uncheckedSize", QC.forAll gen uncheckedSize) :
   ("inBounds", QC.forAll gen inBounds) :
   ("sizeOffset", QC.forAll gen sizeOffset) :
   ("uncheckedSizeOffset", QC.forAll gen uncheckedSizeOffset) :
   ("uncheckedOffset", QC.forAll gen uncheckedOffset) :
   ("lengthIndices", QC.forAll gen lengthIndices) :
   ("indexOffsets", QC.forAll gen indexOffsets) :
   ("invIndices", QC.forAll gen invIndices) :
   ("uncheckedInvIndices", QC.forAll gen uncheckedInvIndices) :
   []
