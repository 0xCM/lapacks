{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Numeric.LAPACK.Matrix.Shape.Private where

import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Extent.Private (Extent)
import Numeric.LAPACK.Wrapper (Flip(Flip, getFlip))

import qualified Type.Data.Num.Unary.Literal as TypeNum
import qualified Type.Data.Num.Unary.Proof as Proof
import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num.Unary (unary, (:+:))
import Type.Data.Num (integralFromProxy)
import Type.Base.Proxy (Proxy(Proxy))

import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Shape (triangleSize, triangleRoot)

import Control.DeepSeq (NFData, rnf)
import Control.Applicative (Const(Const, getConst))

import Data.Functor.Identity (Identity(Identity), runIdentity)
import Data.List (tails)
import Data.Tuple.HT (mapSnd, swap, double)
import Data.Bool.HT (if')


data Order = RowMajor | ColumnMajor
   deriving (Eq, Show)

instance NFData Order where
   rnf RowMajor = ()
   rnf ColumnMajor = ()

flipOrder :: Order -> Order
flipOrder RowMajor = ColumnMajor
flipOrder ColumnMajor = RowMajor

transposeFromOrder :: Order -> Char
transposeFromOrder RowMajor = 'T'
transposeFromOrder ColumnMajor = 'N'

swapOnRowMajor :: Order -> (a,a) -> (a,a)
swapOnRowMajor order =
   case order of
      RowMajor -> swap
      ColumnMajor -> id

sideSwapFromOrder :: Order -> (a,a) -> (Char, (a,a))
sideSwapFromOrder order (m0,n0) =
   let ((side,m), (_,n)) = swapOnRowMajor order (('L', m0), ('R', n0))
   in (side,(m,n))


type family HeightOf shape
type family WidthOf shape


data Full vert horiz height width =
   Full {
      fullOrder :: Order,
      fullExtent :: Extent vert horiz height width
   } deriving (Eq, Show)

type instance HeightOf (Full vert horiz height width) = height
type instance WidthOf (Full vert horiz height width) = width

instance
   (Extent.C vert, Extent.C horiz, NFData height, NFData width) =>
       NFData (Full vert horiz height width) where
   rnf (Full order extent) = rnf (order, extent)

instance
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
      Shape.C (Full vert horiz height width) where

   size (Full _ extent) = Shape.size (Extent.dimensions extent)
   uncheckedSize (Full _ extent) =
      Shape.uncheckedSize (Extent.dimensions extent)

instance
   (Extent.C vert, Extent.C horiz, Shape.Indexed height, Shape.Indexed width) =>
      Shape.Indexed (Full vert horiz height width) where

   type Index (Full vert horiz height width) =
            (Shape.Index height, Shape.Index width)
   indices (Full order extent) = fullIndices order extent

   offset (Full RowMajor extent) =
      Shape.offset (Extent.dimensions extent)
   offset (Full ColumnMajor extent) =
      Shape.offset (swap $ Extent.dimensions extent) . swap
   uncheckedOffset (Full RowMajor extent) =
      Shape.uncheckedOffset (Extent.dimensions extent)
   uncheckedOffset (Full ColumnMajor extent) =
      Shape.uncheckedOffset (swap $ Extent.dimensions extent) . swap

   sizeOffset (Full RowMajor extent) =
      Shape.sizeOffset (Extent.dimensions extent)
   sizeOffset (Full ColumnMajor extent) =
      mapSnd (.swap) $ Shape.sizeOffset (swap $ Extent.dimensions extent)
   uncheckedSizeOffset (Full RowMajor extent) =
      Shape.uncheckedSizeOffset (Extent.dimensions extent)
   uncheckedSizeOffset (Full ColumnMajor extent) =
      mapSnd (.swap) $
      Shape.uncheckedSizeOffset (swap $ Extent.dimensions extent)

   inBounds (Full _ extent) = Shape.inBounds (Extent.dimensions extent)

instance
   (Extent.C vert, Extent.C horiz,
    Shape.InvIndexed height, Shape.InvIndexed width) =>
      Shape.InvIndexed (Full vert horiz height width) where

   indexFromOffset (Full order extent) = fullIndexFromOffset order extent


transpose ::
   (Extent.C vert, Extent.C horiz) =>
   Full vert horiz height width -> Full horiz vert width height
transpose (Full order extent) = Full (flipOrder order) (Extent.transpose extent)

dimensions ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
   Full vert horiz height width -> (Int, Int)
dimensions (Full order extent) =
   swapOnRowMajor order
      (Shape.size $ Extent.height extent,
       Shape.size $ Extent.width extent)

fullHeight ::
   (Extent.C vert, Extent.C horiz) => Full vert horiz height width -> height
fullHeight = Extent.height . fullExtent

fullWidth ::
   (Extent.C vert, Extent.C horiz) => Full vert horiz height width -> width
fullWidth = Extent.width . fullExtent


fullIndices ::
   (Extent.C vert, Extent.C horiz, Shape.Indexed a, Shape.Indexed b) =>
   Order -> Extent vert horiz a b -> [(Shape.Index a, Shape.Index b)]
fullIndices order extent =
   case order of
      RowMajor -> Shape.indices $ Extent.dimensions extent
      ColumnMajor -> map swap $ Shape.indices $ swap $ Extent.dimensions extent

fullIndexFromOffset ::
   (Extent.C vert, Extent.C horiz, Shape.InvIndexed a, Shape.InvIndexed b) =>
   Order -> Extent vert horiz a b -> Int ->
   (Shape.Index a, Shape.Index b)
fullIndexFromOffset order extent =
   case order of
      RowMajor ->
         Shape.indexFromOffset (Extent.dimensions extent)
      ColumnMajor ->
         swap . Shape.indexFromOffset (swap $ Extent.dimensions extent)


type General height width = Full Extent.Big Extent.Big height width
type Tall height width = Full Extent.Big Extent.Small height width
type Wide height width = Full Extent.Small Extent.Big height width
type Square size = Full Extent.Small Extent.Small size size


fullMapExtent ::
   Extent.Map vertA horizA vertB horizB height width ->
   Full vertA horizA height width ->
   Full vertB horizB height width
fullMapExtent f (Full order extent) = Full order $ Extent.apply f extent

general :: Order -> height -> width -> General height width
general order height width = Full order $ Extent.general height width

tall ::
   (Shape.C height, Shape.C width) =>
   Order -> height -> width -> Tall height width
tall order height width =
   if Shape.size height >= Shape.size width
      then Full order $ Extent.tall height width
      else error "MatrixShape.tall: height smaller than width"

wide ::
   (Shape.C height, Shape.C width) =>
   Order -> height -> width -> Wide height width
wide order height width =
   if Shape.size height <= Shape.size width
      then Full order $ Extent.wide height width
      else error "MatrixShape.wide: width smaller than height"

square :: Order -> sh -> Square sh
square order sh = Full order $ Extent.square sh


caseTallWide ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
   Full vert horiz height width ->
   Either (Tall height width) (Wide height width)
caseTallWide (Full order extent) =
   either (Left . Full order) (Right . Full order) $
   Extent.caseTallWide (\h w -> Shape.size h >= Shape.size w) extent


data Split lower vert horiz height width =
   Split {
      splitLower :: lower,
      splitOrder :: Order,
      splitExtent :: Extent vert horiz height width
   } deriving (Eq, Show)

splitHeight ::
   (Extent.C vert, Extent.C horiz) =>
   Split lower vert horiz height width -> height
splitHeight = Extent.height . splitExtent

splitWidth ::
   (Extent.C vert, Extent.C horiz) =>
   Split lower vert horiz height width -> width
splitWidth = Extent.width . splitExtent

splitMapExtent ::
   Extent.Map vertA horizA vertB horizB height width ->
   Split lower vertA horizA height width ->
   Split lower vertB horizB height width
splitMapExtent f (Split lowerPart order extent) =
   Split lowerPart order $ Extent.apply f extent


caseTallWideSplit ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
   Split lower vert horiz height width ->
   Either
      (Split lower Extent.Big Extent.Small height width)
      (Split lower Extent.Small Extent.Big height width)
caseTallWideSplit (Split lowerPart order extent) =
   either (Left . Split lowerPart order) (Right . Split lowerPart order) $
   Extent.caseTallWide (\h w -> Shape.size h >= Shape.size w) extent


type instance HeightOf (Split lower vert horiz height width) = height
type instance WidthOf (Split lower vert horiz height width) = width

data Reflector = Reflector deriving (Eq, Show)
data Triangle = Triangle deriving (Eq, Show)

instance NFData Reflector where rnf Reflector = ()
instance NFData Triangle where rnf Triangle = ()

splitPart ::
   (Extent.C vert, Extent.C horiz,
    Shape.Indexed height, Shape.Indexed width) =>
   Split lower vert horiz height width ->
   (Shape.Index height, Shape.Index width) -> Either lower Triangle
splitPart (Split lowerPart _ extent) (r,c) =
   if Shape.offset (Extent.height extent) r >
         Shape.offset (Extent.width extent) c
     then Left lowerPart
     else Right Triangle

instance
   (NFData lower, Extent.C vert, Extent.C horiz, NFData height, NFData width) =>
      NFData (Split lower vert horiz height width) where
   rnf (Split lowerPart order extent) = rnf (lowerPart, order, extent)

instance
   (Eq lower, Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
      Shape.C (Split lower vert horiz height width) where

   size (Split _ _ extent) = Shape.size (Extent.dimensions extent)
   uncheckedSize (Split _ _ extent) =
      Shape.uncheckedSize (Extent.dimensions extent)

instance
   (Eq lower, Extent.C vert, Extent.C horiz, Shape.Indexed height, Shape.Indexed width) =>
      Shape.Indexed (Split lower vert horiz height width) where

   type Index (Split lower vert horiz height width) =
            (Either lower Triangle,
             (Shape.Index height, Shape.Index width))

   indices sh@(Split _ order extent) =
      map (\ix -> (splitPart sh ix, ix)) $ fullIndices order extent

   offset sh@(Split _ order extent) (part,ix) =
      if part == splitPart sh ix
        then
            case order of
               RowMajor -> Shape.offset (Extent.dimensions extent) ix
               ColumnMajor ->
                  Shape.offset (swap $ Extent.dimensions extent) (swap ix)
        else error "Shape.Split.offset: wrong matrix part"
   uncheckedOffset (Split _ RowMajor extent) =
      Shape.uncheckedOffset (Extent.dimensions extent) . snd
   uncheckedOffset (Split _ ColumnMajor extent) =
      Shape.uncheckedOffset (swap $ Extent.dimensions extent) . swap . snd

   sizeOffset sh@(Split _ order extent) =
      let check (part,ix) a =
            if part == splitPart sh ix
              then a
              else error "Shape.Split.sizeOffset: wrong matrix part"
      in case order of
            RowMajor ->
               mapSnd (\getOffset (part,ix) -> check (part,ix) $ getOffset ix) $
               Shape.sizeOffset (Extent.dimensions extent)
            ColumnMajor ->
               mapSnd
                  (\getOffset (part,ix) ->
                     check (part,ix) $ getOffset $ swap ix) $
               Shape.sizeOffset (swap $ Extent.dimensions extent)
   uncheckedSizeOffset (Split _ RowMajor extent) =
      mapSnd (.snd) $ Shape.uncheckedSizeOffset (Extent.dimensions extent)
   uncheckedSizeOffset (Split _ ColumnMajor extent) =
      mapSnd (.swap.snd) $
      Shape.uncheckedSizeOffset (swap $ Extent.dimensions extent)

   inBounds sh@(Split _ _ extent) (part,ix) =
      Shape.inBounds (Extent.dimensions extent) ix
      &&
      part == splitPart sh ix

instance
   (Eq lower, Extent.C vert, Extent.C horiz,
    Shape.InvIndexed height, Shape.InvIndexed width) =>
      Shape.InvIndexed (Split lower vert horiz height width) where

   indexFromOffset sh@(Split _ order extent) k =
      let ix = fullIndexFromOffset order extent k
      in (splitPart sh ix, ix)


{- |
Store the upper triangular half of a real symmetric or complex Hermitian matrix.
-}
data Hermitian size =
   Hermitian {
      hermitianOrder :: Order,
      hermitianSize :: size
   } deriving (Eq, Show)

type instance HeightOf (Hermitian size) = size
type instance WidthOf (Hermitian size) = size

uploFromOrder :: Order -> Char
uploFromOrder RowMajor = 'L'
uploFromOrder ColumnMajor = 'U'

instance (NFData size) => NFData (Hermitian size) where
   rnf (Hermitian order size) = rnf (order, size)

instance (Shape.C size) => Shape.C (Hermitian size) where
   size (Hermitian _ size) = triangleSize $ Shape.size size
   uncheckedSize (Hermitian _ size) = triangleSize $ Shape.uncheckedSize size

instance (Shape.Indexed size) => Shape.Indexed (Hermitian size) where
   type Index (Hermitian size) = (Shape.Index size, Shape.Index size)

   indices (Hermitian order size) = triangleIndices order size

   offset (Hermitian order size) = triangleOffset order size
   uncheckedOffset (Hermitian order size) =
      triangleUncheckedOffset order size

   sizeOffset (Hermitian order size) = triangleSizeOffset order size
   uncheckedSizeOffset (Hermitian order size) =
      triangleUncheckedSizeOffset order size

   inBounds (Hermitian _ size) ix@(r,c) =
      Shape.inBounds (size,size) ix
      &&
      Shape.offset size r <= Shape.offset size c

instance (Shape.InvIndexed size) => Shape.InvIndexed (Hermitian size) where
   indexFromOffset (Hermitian order size) k =
      triangleIndexFromOffset order size k



data Triangular lo diag up size =
   Triangular {
      triangularDiag :: diag,
      triangularUplo :: (lo,up),
      triangularOrder :: Order,
      triangularSize :: size
   } deriving (Eq, Show)

type instance HeightOf (Triangular lo diag up size) = size
type instance WidthOf (Triangular lo diag up size) = size


data Unit = Unit deriving (Eq, Show)
data NonUnit = NonUnit deriving (Eq, Show)

class TriDiag diag where switchTriDiag :: f Unit -> f NonUnit -> f diag
instance TriDiag Unit where switchTriDiag f _ = f
instance TriDiag NonUnit where switchTriDiag _ f = f

autoDiag :: TriDiag diag => diag
autoDiag = runIdentity $ switchTriDiag (Identity Unit) (Identity NonUnit)

caseTriDiag :: TriDiag diag => diag -> a -> a -> a
caseTriDiag diag unit nonUnit =
   getConstAs diag $ switchTriDiag (Const unit) (Const nonUnit)

charFromTriDiag :: TriDiag diag => diag -> Char
charFromTriDiag diag = caseTriDiag diag 'U' 'N'


relaxUnitDiagonal ::
   (TriDiag diag) => Triangular lo Unit up sh -> Triangular lo diag up sh
relaxUnitDiagonal shape = shape{triangularDiag = autoDiag}

strictNonUnitDiagonal ::
   (TriDiag diag) => Triangular lo diag up sh -> Triangular lo NonUnit up sh
strictNonUnitDiagonal shape = shape{triangularDiag = NonUnit}


data Empty = Empty deriving (Eq, Show)
data Filled = Filled deriving (Eq, Show)

lower :: (Filled,Empty)
lower = (Filled,Empty)

upper :: (Empty,Filled)
upper = (Empty,Filled)

type Identity = Triangular Empty Unit Empty
type Diagonal = Triangular Empty NonUnit Empty
type LowerTriangular diag = Triangular Filled diag Empty
type UpperTriangular diag = Triangular Empty diag Filled
type FlexSymmetric diag = Triangular Filled diag Filled
type Symmetric = FlexSymmetric NonUnit

triangularTranspose ::
   (Content lo, Content up) =>
   Triangular lo diag up sh -> Triangular up diag lo sh
triangularTranspose (Triangular diag uplo order size) =
   Triangular diag
      (swap uplo)
      (caseDiagUpLoSym uplo flipOrder flipOrder flipOrder id order)
      size


class Content c where switchContent :: f Empty -> f Filled -> f c
instance Content Empty where switchContent f _ = f
instance Content Filled where switchContent _ f = f


type UpLo lo up = (UpLoC lo up, UpLoC up lo)

class (DiagUpLoC lo up, UpLoSymC lo up) => UpLoC lo up where
   switchUpLo :: f Empty Filled -> f Filled Empty -> f lo up

instance UpLoC Empty  Filled where switchUpLo f _ = f
instance UpLoC Filled Empty  where switchUpLo _ f = f


type DiagUpLo lo up = (DiagUpLoC lo up, DiagUpLoC up lo)

class (Content lo, Content up) => DiagUpLoC lo up where
   switchDiagUpLo ::
      f Empty Empty -> f Empty Filled -> f Filled Empty -> f lo up

instance DiagUpLoC Empty  Empty  where switchDiagUpLo f _ _ = f
instance DiagUpLoC Empty  Filled where switchDiagUpLo _ f _ = f
instance DiagUpLoC Filled Empty  where switchDiagUpLo _ _ f = f


type UpLoSym lo up = (UpLoSymC lo up, UpLoSymC up lo)

class (Content lo, Content up) => UpLoSymC lo up where
   switchUpLoSym ::
      f Empty Filled -> f Filled Empty -> f Filled Filled -> f lo up

instance UpLoSymC Empty  Filled where switchUpLoSym f _ _ = f
instance UpLoSymC Filled Empty  where switchUpLoSym _ f _ = f
instance UpLoSymC Filled Filled where switchUpLoSym _ _ f = f


switchDiagUpLoSym ::
   (Content lo, Content up) =>
   f Empty Empty -> f Empty Filled -> f Filled Empty -> f Filled Filled ->
   f lo up
switchDiagUpLoSym fDiag fUpper fLower fSymm =
   getFlip $
   switchContent
      (Flip $ switchContent fDiag fUpper)
      (Flip $ switchContent fLower fSymm)

autoContent :: Content c => c
autoContent = runIdentity $ switchContent (Identity Empty) (Identity Filled)

autoUplo :: (Content lo, Content up) => (lo,up)
autoUplo = (autoContent,autoContent)

uploOrder :: (Content lo, Content up) => (lo,up) -> Order -> Order
uploOrder (_loc,upc) = caseContent upc flipOrder id

getConstAs :: c -> Const a c -> a
getConstAs _ = getConst

caseContent :: Content c => c -> a -> a -> a
caseContent c lo up =
   getConstAs c $ switchContent (Const lo) (Const up)

caseLoUp :: UpLo lo up => (lo,up) -> a -> a -> a
caseLoUp (_loc,upc) = caseContent upc

caseDiagUpLoSym :: (Content lo, Content up) => (lo,up) -> a -> a -> a -> a -> a
caseDiagUpLoSym (loc,upc) diag up lo symm =
   caseContent loc
      (caseContent upc diag up)
      (caseContent upc lo symm)


newtype Const2 a lo up = Const2 {getConst2 :: a}

getContentConst2 :: (lo,up) -> Const2 a lo up -> a
getContentConst2 _ = getConst2

caseUpLoSym :: (UpLoSym lo up) => (lo,up) -> a -> a -> a -> a
caseUpLoSym c lo up sym =
   getContentConst2 c $ switchUpLoSym (Const2 lo) (Const2 up) (Const2 sym)


instance
   (Content lo, TriDiag diag, Content up, NFData size) =>
      NFData (Triangular lo diag up size) where
   rnf (Triangular diag (loc,upc) order size) =
      rnf
         (flip getFlip diag $
            switchTriDiag (Flip $ \Unit -> ()) (Flip $ \NonUnit -> ()),
          let rnfContent c =
               flip getFlip c $
               switchContent
                  (Flip $ \Empty -> ())
                  (Flip $ \Filled -> ())
          in (rnfContent loc, rnfContent upc),
          order, size)

instance
   (Content lo, TriDiag diag, Content up, Shape.C size) =>
      Shape.C (Triangular lo diag up size) where

   size (Triangular _diag uplo _ size) =
      let n = Shape.size size
      in caseDiagUpLoSym uplo n
            (triangleSize n)
            (triangleSize n)
            (triangleSize n)
   uncheckedSize (Triangular _diag uplo _ size) =
      let n = Shape.uncheckedSize size
      in caseDiagUpLoSym uplo n
            (triangleSize n)
            (triangleSize n)
            (triangleSize n)

instance
   (Content lo, TriDiag diag, Content up, Shape.Indexed size) =>
      Shape.Indexed (Triangular lo diag up size) where
   type Index (Triangular lo diag up size) =
         (Shape.Index size, Shape.Index size)

   indices (Triangular _diag uplo order size) =
      caseDiagUpLoSym uplo
         (map double $ Shape.indices size)
         (triangleIndices order size)
         (map swap $ triangleIndices (flipOrder order) size)
         (triangleIndices order size)

   offset (Triangular _diag uplo order size) =
      caseDiagUpLoSym uplo
         (Shape.offset size . snd)
         (triangleOffset order size)
         (triangleOffset (flipOrder order) size . swap)
         (triangleOffset order size)
   uncheckedOffset (Triangular _diag uplo order size) =
      caseDiagUpLoSym uplo
         (Shape.offset size . snd)
         (triangleUncheckedOffset order size)
         (triangleUncheckedOffset (flipOrder order) size . swap)
         (triangleUncheckedOffset order size)

   sizeOffset (Triangular _diag uplo order size) =
      caseDiagUpLoSym uplo
         (mapSnd (.snd) $ Shape.sizeOffset size)
         (triangleSizeOffset order size)
         (mapSnd (.swap) $ triangleSizeOffset (flipOrder order) size)
         (triangleSizeOffset order size)
   uncheckedSizeOffset (Triangular _diag uplo order size) =
      caseDiagUpLoSym uplo
         (mapSnd (.snd) $ Shape.uncheckedSizeOffset size)
         (triangleUncheckedSizeOffset order size)
         (mapSnd (.swap) $ triangleUncheckedSizeOffset (flipOrder order) size)
         (triangleUncheckedSizeOffset order size)

   inBounds (Triangular _diag uplo _ size) ix@(r,c) =
      Shape.inBounds (size,size) ix
      &&
      caseDiagUpLoSym uplo
         (Shape.offset size r == Shape.offset size c)
         (Shape.offset size r <= Shape.offset size c)
         (Shape.offset size r >= Shape.offset size c)
         (Shape.offset size r <= Shape.offset size c)

instance
   (Content lo, TriDiag diag, Content up, Shape.InvIndexed size) =>
      Shape.InvIndexed (Triangular lo diag up size) where

   indexFromOffset (Triangular _diag uplo order size) k =
      caseDiagUpLoSym uplo
         (double $ Shape.indexFromOffset size k)
         (triangleIndexFromOffset order size k)
         (swap $ triangleIndexFromOffset (flipOrder order) size k)
         (triangleIndexFromOffset order size k)


triangleRootDouble :: Int -> Double
triangleRootDouble = triangleRoot . fromIntegral

triangleExtent :: String -> Int -> Int
triangleExtent name size =
   let n = round $ triangleRootDouble size
   in if size == triangleSize n
        then n
        else error (name ++ ": no triangular number of elements")

triangleIndices ::
   (Shape.Indexed sh) => Order -> sh -> [(Shape.Index sh, Shape.Index sh)]
triangleIndices RowMajor = Shape.indices . Shape.upperTriangular
triangleIndices ColumnMajor = map swap . Shape.indices . Shape.lowerTriangular

triangleOffset ::
   (Shape.Indexed sh) => Order -> sh -> (Shape.Index sh, Shape.Index sh) -> Int
triangleOffset order size =
   case order of
      RowMajor -> Shape.offset (Shape.upperTriangular size)
      ColumnMajor -> Shape.offset (Shape.lowerTriangular size) . swap

triangleUncheckedOffset ::
   (Shape.Indexed sh) => Order -> sh -> (Shape.Index sh, Shape.Index sh) -> Int
triangleUncheckedOffset order size =
   case order of
      RowMajor -> Shape.uncheckedOffset (Shape.upperTriangular size)
      ColumnMajor -> Shape.uncheckedOffset (Shape.lowerTriangular size) . swap

triangleSizeOffset ::
   (Shape.Indexed sh) =>
   Order -> sh -> (Int, (Shape.Index sh, Shape.Index sh) -> Int)
triangleSizeOffset order size =
   case order of
      RowMajor -> Shape.sizeOffset (Shape.upperTriangular size)
      ColumnMajor ->
         mapSnd (.swap) $ Shape.sizeOffset (Shape.lowerTriangular size)

triangleUncheckedSizeOffset ::
   (Shape.Indexed sh) =>
   Order -> sh -> (Int, (Shape.Index sh, Shape.Index sh) -> Int)
triangleUncheckedSizeOffset order size =
   case order of
      RowMajor -> Shape.uncheckedSizeOffset (Shape.upperTriangular size)
      ColumnMajor ->
         mapSnd (.swap) $ Shape.uncheckedSizeOffset (Shape.lowerTriangular size)

triangleIndexFromOffset ::
   (Shape.InvIndexed sh) =>
   Order -> sh -> Int -> (Shape.Index sh, Shape.Index sh)
triangleIndexFromOffset order size =
   case order of
      RowMajor -> Shape.indexFromOffset (Shape.upperTriangular size)
      ColumnMajor -> swap . Shape.indexFromOffset (Shape.lowerTriangular size)


type UnaryProxy a = Proxy (Unary.Un a)

data Banded sub super vert horiz height width =
   Banded {
      bandedOffDiagonals :: (UnaryProxy sub, UnaryProxy super),
      bandedOrder :: Order,
      bandedExtent :: Extent vert horiz height width
   } deriving (Eq, Show)

type BandedGeneral sub super = Banded sub super Extent.Big Extent.Big
type BandedSquare sub super size =
      Banded sub super Extent.Small Extent.Small size size

type BandedLowerTriangular sub size = BandedSquare sub TypeNum.U0 size
type BandedUpperTriangular super size = BandedSquare TypeNum.U0 super size

type BandedDiagonal size = BandedSquare TypeNum.U0 TypeNum.U0 size


bandedHeight ::
   (Extent.C vert, Extent.C horiz) =>
   Banded sub super vert horiz height width -> height
bandedHeight = Extent.height . bandedExtent

bandedWidth ::
   (Extent.C vert, Extent.C horiz) =>
   Banded sub super vert horiz height width -> width
bandedWidth = Extent.width . bandedExtent

bandedMapExtent ::
   Extent.Map vertA horizA vertB horizB height width ->
   Banded sub super vertA horizA height width ->
   Banded sub super vertB horizB height width
bandedMapExtent f (Banded numDiag order extent) =
   Banded numDiag order $ Extent.apply f extent

type instance HeightOf (Banded sub super vert horiz height width) = height
type instance WidthOf (Banded sub super vert horiz height width) = width

bandedBreadth ::
   (Unary.Natural sub, Unary.Natural super) =>
   (UnaryProxy sub, UnaryProxy super) -> Int
bandedBreadth (sub,super) =
   integralFromProxy sub + 1 + integralFromProxy super

numOffDiagonals ::
   (Unary.Natural sub, Unary.Natural super) =>
   Order -> (UnaryProxy sub, UnaryProxy super) -> (Int,Int)
numOffDiagonals order (sub,super) =
   swapOnRowMajor order (integralFromProxy sub, integralFromProxy super)

natFromProxy :: (Unary.Natural n) => UnaryProxy n -> Proof.Nat n
natFromProxy Proxy = Proof.Nat

addOffDiagonals ::
   (Unary.Natural subA, Unary.Natural superA,
    Unary.Natural subB, Unary.Natural superB,
    (subA :+: subB) ~ subC,
    (superA :+: superB) ~ superC) =>
   (UnaryProxy subA, UnaryProxy superA) ->
   (UnaryProxy subB, UnaryProxy superB) ->
   ((Proof.Nat subC, Proof.Nat superC),
    (UnaryProxy subC, UnaryProxy superC))
addOffDiagonals (subA,superA) (subB,superB) =
   ((Proof.addNat (natFromProxy subA) (natFromProxy subB),
     Proof.addNat (natFromProxy superA) (natFromProxy superB)),
    (Proxy,Proxy))

bandedTranspose ::
   (Extent.C vert, Extent.C horiz) =>
   Banded sub super vert horiz height width ->
   Banded super sub horiz vert width height
bandedTranspose (Banded (sub,super) order extent) =
   Banded (super,sub) (flipOrder order) (Extent.transpose extent)


bandedGeneral ::
   (UnaryProxy sub, UnaryProxy super) -> Order -> height -> width ->
   Banded sub super Extent.Big Extent.Big height width
bandedGeneral offDiag order height width =
   Banded offDiag order (Extent.general height width)

bandedSquare ::
   (UnaryProxy sub, UnaryProxy super) -> Order -> size ->
   Banded sub super Extent.Small Extent.Small size size
bandedSquare offDiag order = Banded offDiag order . Extent.square


data BandedIndex row column =
     InsideBox row column
   | VertOutsideBox Int column
   | HorizOutsideBox row Int
   deriving (Eq, Show)

instance
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz, NFData height, NFData width) =>
      NFData (Banded sub super vert horiz height width) where
   rnf (Banded (Proxy,Proxy) order extent) = rnf (order, extent)

instance
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
      Shape.C (Banded sub super vert horiz height width) where

   size (Banded offDiag order extent) =
      bandedBreadth offDiag *
      case order of
         RowMajor -> Shape.size (Extent.height extent)
         ColumnMajor -> Shape.size (Extent.width extent)
   uncheckedSize (Banded offDiag order extent) =
      bandedBreadth offDiag *
      case order of
         RowMajor -> Shape.uncheckedSize (Extent.height extent)
         ColumnMajor -> Shape.uncheckedSize (Extent.width extent)

instance
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz, Shape.Indexed height, Shape.Indexed width) =>
      Shape.Indexed (Banded sub super vert horiz height width) where

   type Index (Banded sub super vert horiz height width) =
            BandedIndex (Shape.Index height) (Shape.Index width)
   indices (Banded (sub,super) order extent) =
      let (height,width) = Extent.dimensions extent
      in case order of
            RowMajor ->
               map (\(r,c) -> either (HorizOutsideBox r) (InsideBox r) c) $
               bandedIndicesRowMajor (sub,super) (height,width)
            ColumnMajor ->
               map (\(c,r) ->
                     either (flip VertOutsideBox c) (flip InsideBox c) r) $
               bandedIndicesRowMajor (super,sub) (width,height)

   offset shape ix =
      if Shape.inBounds shape ix
         then Shape.uncheckedOffset shape ix
         else error "Banded.offset: index outside band"

   uncheckedOffset (Banded (sub,super) order extent) ix =
      let (height,width) = Extent.dimensions extent
          kl = integralFromProxy sub
          ku = integralFromProxy super
      in bandedOffset (kl,ku) order (height,width) ix

   inBounds (Banded (sub,super) order extent) ix =
      let (height,width) = Extent.dimensions extent
          kl = integralFromProxy sub
          ku = integralFromProxy super
          insideBand r c = Shape.inBounds (Shape.Range (-kl) ku) (c-r)
      in case (order,ix) of
            (_, InsideBox r c) ->
               Shape.inBounds (height,width) (r,c)
               &&
               insideBand (Shape.offset height r) (Shape.offset width c)
            (RowMajor, HorizOutsideBox r c) ->
               Shape.inBounds height r
               &&
               insideBand (Shape.offset height r) (outsideOffset width c)
            (ColumnMajor, VertOutsideBox r c) ->
               Shape.inBounds width c
               &&
               insideBand (outsideOffset height r) (Shape.offset width c)
            _ -> False

instance
   (Unary.Natural sub, Unary.Natural super, Extent.C vert, Extent.C horiz,
    Shape.InvIndexed height, Shape.InvIndexed width) =>
      Shape.InvIndexed (Banded sub super vert horiz height width) where

   indexFromOffset (Banded (sub,super) order extent) j =
      bandedIndexFromOffset
         Shape.indexFromOffset Shape.indexFromOffset
         (integralFromProxy sub, integralFromProxy super) order
         (Extent.dimensions extent) j

   uncheckedIndexFromOffset (Banded (sub,super) order extent) j =
      bandedIndexFromOffset
         Shape.uncheckedIndexFromOffset Shape.uncheckedIndexFromOffset
         (integralFromProxy sub, integralFromProxy super) order
         (Extent.dimensions extent) j

outsideOffset :: Shape.C sh => sh -> Int -> Int
outsideOffset size k = if k<0 then k else Shape.size size + k

bandedOffset ::
   (Shape.Indexed height, Shape.Indexed width) =>
   (Int, Int) -> Order -> (height, width) ->
   BandedIndex (Shape.Index height) (Shape.Index width) -> Int
bandedOffset (kl,ku) order (height,width) ix =
   let k = kl+ku
   in case ix of
         InsideBox r c ->
            let i = Shape.uncheckedOffset height r
                j = Shape.uncheckedOffset width c
            in case order of
                  RowMajor -> k*i + kl+j
                  ColumnMajor -> k*j + ku+i
         VertOutsideBox r c ->
            let i = outsideOffset height r
                j = Shape.uncheckedOffset width c
            in  k*j + ku+i
         HorizOutsideBox r c ->
            let i = Shape.uncheckedOffset height r
                j = outsideOffset width c
            in  k*i + kl+j

bandedIndicesRowMajor ::
   (Unary.Natural sub, Unary.Natural super,
    Shape.Indexed height, Shape.Indexed width) =>
   (UnaryProxy sub, UnaryProxy super) ->
   (height, width) ->
   [(Shape.Index height, Either Int (Shape.Index width))]
bandedIndicesRowMajor (sub,super) (height,width) =
   let kl = integralFromProxy sub
       ku = integralFromProxy super
   in concat $
      zipWith (\r -> map ((,) r)) (Shape.indices height) $
      map (take (kl+1+ku)) $ tails $
         (map Left $ take kl $ iterate (1+) (-kl)) ++
         (map Right $ Shape.indices width) ++
         (map Left $ iterate (1+) 0)

bandedIndexFromOffset ::
   (Shape.C height, Shape.C width) =>
   (height -> Int -> row) ->
   (width -> Int -> column) ->
   (Int,Int) -> Order -> (height,width) -> Int -> BandedIndex row column
bandedIndexFromOffset
      rowFromOffset columnFromOffset (kl,ku) order (height,width) j =
   case order of
      RowMajor ->
         let n = Shape.size width
             (rb,cb) = divMod j (kl+1+ku)
             r = rowFromOffset height rb
             ci = rb+cb-kl
         in if' (ci<0) (HorizOutsideBox r ci) $
            if' (ci>=n) (HorizOutsideBox r (ci-n)) $
            InsideBox r (columnFromOffset width ci)
      ColumnMajor ->
         let m = Shape.size height
             (cb,rb) = divMod j (kl+1+ku)
             c = columnFromOffset width cb
             ri = rb+cb-ku
         in if' (ri<0) (VertOutsideBox ri c) $
            if' (ri>=m) (VertOutsideBox (ri-m) c) $
            InsideBox (rowFromOffset height ri) c


data BandedHermitian off size =
   BandedHermitian {
      bandedHermitianOffDiagonals :: UnaryProxy off,
      bandedHermitianOrder :: Order,
      bandedHermitianSize :: size
   } deriving (Eq, Show)

type instance HeightOf (BandedHermitian off size) = size
type instance WidthOf (BandedHermitian off size) = size

instance (Unary.Natural off, NFData size) =>
      NFData (BandedHermitian off size) where
   rnf (BandedHermitian Proxy order size) = rnf (order, size)

instance (Unary.Natural off, Shape.C size) =>
      Shape.C (BandedHermitian off size) where
   size (BandedHermitian offDiag _order size) =
      (1 + integralFromProxy offDiag) * Shape.size size
   uncheckedSize (BandedHermitian offDiag _order size) =
      (1 + integralFromProxy offDiag) * Shape.uncheckedSize size

instance (Unary.Natural off, Shape.Indexed size) =>
      Shape.Indexed (BandedHermitian off size) where
   type Index (BandedHermitian off size) =
            BandedIndex (Shape.Index size) (Shape.Index size)
   indices (BandedHermitian offDiag order size) =
      case order of
         RowMajor ->
            map (\(r,c) -> either (HorizOutsideBox r) (InsideBox r) c) $
            bandedIndicesRowMajor (unary TypeNum.u0, offDiag) (size,size)
         ColumnMajor ->
            map (\(c,r) ->
                  either (flip VertOutsideBox c) (flip InsideBox c) r) $
            bandedIndicesRowMajor (offDiag, unary TypeNum.u0) (size,size)

   offset shape ix =
      if Shape.inBounds shape ix
         then Shape.uncheckedOffset shape ix
         else error "BandedHermitian.offset: index outside band"

   uncheckedOffset (BandedHermitian offDiag order size) ix =
      let k = integralFromProxy offDiag
      in bandedOffset (0,k) order (size,size) ix

   inBounds (BandedHermitian offDiag order size) ix =
      let ku = integralFromProxy offDiag
          insideBand r c = Shape.inBounds (Shape.Range 0 ku) (c-r)
      in case (order,ix) of
            (_, InsideBox r c) ->
               Shape.inBounds (size,size) (r,c)
               &&
               insideBand (Shape.offset size r) (Shape.offset size c)
            (RowMajor, HorizOutsideBox r c) ->
               Shape.inBounds size r
               &&
               insideBand (Shape.offset size r) (outsideOffset size c)
            (ColumnMajor, VertOutsideBox r c) ->
               Shape.inBounds size c
               &&
               insideBand (outsideOffset size r) (Shape.offset size c)
            _ -> False

instance (Unary.Natural off, Shape.InvIndexed size) =>
      Shape.InvIndexed (BandedHermitian off size) where

   indexFromOffset (BandedHermitian offDiag order size) j =
      bandedHermitianIndexFromOffset
         Shape.indexFromOffset Shape.indexFromOffset
         (integralFromProxy offDiag) order size j

   uncheckedIndexFromOffset (BandedHermitian offDiag order size) j =
      bandedHermitianIndexFromOffset
         Shape.uncheckedIndexFromOffset Shape.uncheckedIndexFromOffset
         (integralFromProxy offDiag) order size j

bandedHermitianIndexFromOffset ::
   (Shape.C sh) =>
   (sh -> Int -> row) ->
   (sh -> Int -> column) ->
   Int -> Order -> sh -> Int -> BandedIndex row column
bandedHermitianIndexFromOffset rowFromOffset columnFromOffset k order size j =
   case order of
      RowMajor ->
         let n = Shape.size size
             (rb,cb) = divMod j (k+1)
             r = rowFromOffset size rb
             ci = rb+cb
         in if ci<n
               then InsideBox r (columnFromOffset size ci)
               else HorizOutsideBox r (ci-n)
      ColumnMajor ->
         let (cb,rb) = divMod j (k+1)
             c = columnFromOffset size cb
             ri = rb+cb-k
         in if ri>=0
               then InsideBox (rowFromOffset size ri) c
               else VertOutsideBox ri c
