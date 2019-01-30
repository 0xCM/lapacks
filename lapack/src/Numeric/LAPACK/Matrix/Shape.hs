module Numeric.LAPACK.Matrix.Shape (
   General,
   Tall,
   Wide,
   Square,
   Full(..), fullHeight, fullWidth,
   Order(..), flipOrder,
   general,
   square,
   wide,
   tall,

   Split,
   SplitGeneral,
   Triangle(..),
   Reflector(..),
   splitGeneral,
   splitFromFull,

   Hermitian(..),
   hermitian,

   Triangular(..),
   Identity,
   Diagonal,
   LowerTriangular,
   UpperTriangular,
   Symmetric,
   diagonal,
   lowerTriangular,
   upperTriangular,
   symmetric,
   autoDiag,
   autoUplo,
   DiagUpLo,
   Unit(Unit),
   NonUnit(NonUnit),

   Banded(..),
   BandedGeneral,
   BandedSquare,
   BandedLowerTriangular,
   BandedUpperTriangular,
   BandedDiagonal,
   BandedIndex(..),
   bandedGeneral,
   bandedSquare,
   bandedFromFull,
   UnaryProxy,
   addOffDiagonals,
   TriDiag,
   switchTriDiag,
   Content,

   BandedHermitian(..),
   bandedHermitian,
   ) where

import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Shape.Private


type SplitGeneral lower height width =
      Split lower Extent.Big Extent.Big height width

splitGeneral ::
   lower -> Order -> height -> width -> SplitGeneral lower height width
splitGeneral lowerPart order height width =
   Split lowerPart order $ Extent.general height width

splitFromFull ::
   lower ->
   Full vert horiz height width ->
   Split lower vert horiz height width
splitFromFull lowerPart (Full order extent) = Split lowerPart order extent


diagonal :: Order -> size -> Triangular Empty NonUnit Empty size
diagonal = Triangular NonUnit autoUplo

lowerTriangular :: Order -> size -> LowerTriangular NonUnit size
lowerTriangular = Triangular NonUnit autoUplo

upperTriangular :: Order -> size -> UpperTriangular NonUnit size
upperTriangular = Triangular NonUnit autoUplo

symmetric :: Order -> size -> Symmetric size
symmetric = Triangular NonUnit autoUplo

hermitian :: Order -> size -> Hermitian size
hermitian = Hermitian


bandedFromFull ::
   (UnaryProxy sub, UnaryProxy super) ->
   Full vert horiz height width ->
   Banded sub super vert horiz height width
bandedFromFull offDiag (Full order extent) = Banded offDiag order extent


bandedHermitian :: UnaryProxy off -> Order -> size -> BandedHermitian off size
bandedHermitian = BandedHermitian
