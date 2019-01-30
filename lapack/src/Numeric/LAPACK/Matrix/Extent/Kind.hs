{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE GADTs #-}
module Numeric.LAPACK.Matrix.Extent.Kind where


data General height width =
   General {
      generalHeight :: height,
      generalWidth :: width
   } deriving (Eq, Show)

data Tall height width =
   Tall {
      tallHeight :: height,
      tallWidth :: width
   } deriving (Eq, Show)

data Wide height width =
   Wide {
      wideHeight :: height,
      wideWidth :: width
   } deriving (Eq, Show)

data Square height width =
   (height ~ width) =>
   Square {
      squareSize :: height
   }

instance (Eq height, Eq width) => Eq (Square height width) where
   Square a == Square b  =  a==b

instance (Show height, Show width) => Show (Square height width) where
   showsPrec p (Square s) =
      showParen (p>10) (showString "Square " . showsPrec 11 s)
