{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Numeric.LAPACK.Matrix.Banded (
   module Numeric.LAPACK.Matrix.Banded.Basic,
   height, width,

   solve,
   determinant,
   ) where

import Numeric.LAPACK.Matrix.Banded.Basic
import Numeric.LAPACK.Matrix.Banded.Linear

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent

import qualified Data.Array.Comfort.Storable as Array


height ::
   (Extent.C vert, Extent.C horiz) =>
   Banded sub super vert horiz height width a -> height
height = MatrixShape.bandedHeight . Array.shape

width ::
   (Extent.C vert, Extent.C horiz) =>
   Banded sub super vert horiz height width a -> width
width = MatrixShape.bandedWidth . Array.shape
