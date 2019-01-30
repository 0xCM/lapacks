module Numeric.LAPACK.Matrix.Private where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent as Extent
import Numeric.LAPACK.Matrix.Shape.Private (Order, flipOrder)

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.ForeignPtr (ForeignPtr)


type Full vert horiz height width =
         Array (MatrixShape.Full vert horiz height width)

type General height width = Array (MatrixShape.General height width)
type Tall height width = Array (MatrixShape.Tall height width)
type Wide height width = Array (MatrixShape.Wide height width)
type Square sh = Array (MatrixShape.Square sh)


argGeneral ::
   (MatrixShape.Order -> height -> width -> ForeignPtr a -> b) ->
   (General height width a -> b)
argGeneral f (Array (MatrixShape.Full order extent) a) =
   f order (Extent.height extent) (Extent.width extent) a

argSquare ::
   (MatrixShape.Order -> sh -> ForeignPtr a -> b) -> (Square sh a -> b)
argSquare f (Array (MatrixShape.Full order extent) a) =
   f order (Extent.squareSize extent) a


type ZeroInt = Shape.ZeroBased Int

zeroInt :: Int -> ZeroInt
zeroInt = Shape.ZeroBased


mapExtent ::
   (Extent.C vertA, Extent.C horizA) =>
   (Extent.C vertB, Extent.C horizB) =>
   Extent.Map vertA horizA vertB horizB height width ->
   Full vertA horizA height width a -> Full vertB horizB height width a
mapExtent f = Array.mapShape $ MatrixShape.fullMapExtent f

fromFull ::
   (Extent.C vert, Extent.C horiz) =>
   Full vert horiz height width a -> General height width a
fromFull = mapExtent Extent.toGeneral

generalizeTall ::
   (Extent.C vert, Extent.C horiz) =>
   Full vert Extent.Small height width a -> Full vert horiz height width a
generalizeTall = mapExtent Extent.generalizeTall

generalizeWide ::
   (Extent.C vert, Extent.C horiz) =>
   Full Extent.Small horiz height width a -> Full vert horiz height width a
generalizeWide = mapExtent Extent.generalizeWide


height ::
   (Extent.C vert, Extent.C horiz) =>
   Full vert horiz height width a -> height
height = MatrixShape.fullHeight . Array.shape

width ::
   (Extent.C vert, Extent.C horiz) =>
   Full vert horiz height width a -> width
width = MatrixShape.fullWidth . Array.shape


data Transposition = NonTransposed | Transposed
   deriving (Eq, Show, Enum, Bounded)

transposeOrder :: Transposition -> Order -> Order
transposeOrder NonTransposed = id
transposeOrder Transposed = flipOrder

data Conjugation = NonConjugated | Conjugated
   deriving (Eq, Show, Enum, Bounded)

data Inversion = NonInverted | Inverted
   deriving (Eq, Show, Enum, Bounded)

flipInversion :: Inversion -> Inversion
flipInversion NonInverted = Inverted
flipInversion Inverted = NonInverted
