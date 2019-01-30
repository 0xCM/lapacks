module Numeric.LAPACK.Matrix.Extent (
   Extent.C(switchTag),
   Extent.Extent,
   Map,
   Small, Big,
   Extent.height,
   Extent.width,
   Extent.squareSize,
   Extent.dimensions,
   Extent.transpose,
   Extent.fuse,

   Extent.square,

   toGeneral,
   fromSquare,
   fromSquareLiberal,
   generalizeTall,
   generalizeWide,
   ) where

import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Extent.Private (C, Small, Big, Map(Map))


toGeneral ::
   (C vert, C horiz) => Map vert horiz Big Big height width
toGeneral = Map Extent.toGeneral

fromSquare :: (C vert, C horiz) => Map Small Small vert horiz size size
fromSquare = Map Extent.fromSquare

fromSquareLiberal ::
   (C vert, C horiz) => Map Small Small vert horiz height width
fromSquareLiberal = Map Extent.fromSquareLiberal

generalizeTall :: (C vert, C horiz) => Map vert Small vert horiz height width
generalizeTall = Map Extent.generalizeTall

generalizeWide :: (C vert, C horiz) => Map Small horiz vert horiz height width
generalizeWide = Map Extent.generalizeWide
