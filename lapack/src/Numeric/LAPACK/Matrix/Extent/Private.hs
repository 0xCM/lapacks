{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE GADTs #-}
module Numeric.LAPACK.Matrix.Extent.Private where

import qualified Numeric.LAPACK.Matrix.Extent.Kind as EK
import Numeric.LAPACK.Wrapper (Flip(Flip, getFlip))

import Control.DeepSeq (NFData, rnf)

import Data.Maybe.HT (toMaybe)
import Data.Tuple.HT (swap)
import Data.Eq.HT (equating)


data Extent vertical horizontal height width =
   Extent {
      extentDir :: (vertical,horizontal),
      extentDim :: Dimensions vertical horizontal height width
   }

instance
   (C vertical, C horizontal, NFData height, NFData width) =>
      NFData (Extent vertical horizontal height width) where
   rnf =
      getAccessor $
      switchTagPair
         (Accessor $ \(Extent o (EK.Square s)) -> rnf (o,s))
         (Accessor $ \(Extent o (EK.Wide h w)) -> rnf (o,(h,w)))
         (Accessor $ \(Extent o (EK.Tall h w)) -> rnf (o,(h,w)))
         (Accessor $ \(Extent o (EK.General h w)) -> rnf (o,(h,w)))


data Big = Big deriving (Eq,Show)
data Small = Small deriving (Eq,Show)

instance NFData Big where rnf Big = ()
instance NFData Small where rnf Small = ()

type General = Extent Big Big
type Tall = Extent Big Small
type Wide = Extent Small Big
type Square sh = Extent Small Small sh sh


type family Dimensions vertical horizontal :: * -> * -> *

type instance Dimensions Big Big = EK.General
type instance Dimensions Big Small = EK.Tall
type instance Dimensions Small Big = EK.Wide
type instance Dimensions Small Small = EK.Square


general :: height -> width -> General height width
general h w = Extent (Big,Big) $ EK.General h w

tall :: height -> width -> Tall height width
tall h w = Extent (Big,Small) $ EK.Tall h w

wide :: height -> width -> Wide height width
wide h w = Extent (Small,Big) $ EK.Wide h w

square :: sh -> Square sh
square sh = Extent (Small,Small) $ EK.Square sh


newtype Map vertA horizA vertB horizB height width =
   Map {
      apply ::
         Extent vertA horizA height width ->
         Extent vertB horizB height width
   }


class C tag where switchTag :: f Small -> f Big -> f tag
instance C Small where switchTag f _ = f
instance C Big where switchTag _ f = f


switchTagPair ::
   (C vert, C horiz) =>
   f Small Small -> f Small Big -> f Big Small -> f Big Big -> f vert horiz
switchTagPair fSquare fWide fTall fGeneral =
   getFlip $
   switchTag
      (Flip $ switchTag fSquare fWide)
      (Flip $ switchTag fTall fGeneral)


newtype CaseTallWide height width vert horiz =
   CaseTallWide {
      getCaseTallWide ::
         Extent vert horiz height width ->
         Either (Tall height width) (Wide height width)
   }

caseTallWide ::
   (C vert, C horiz) =>
   (height -> width -> Bool) ->
   Extent vert horiz height width ->
   Either (Tall height width) (Wide height width)
caseTallWide ge =
   getCaseTallWide $
   switchTagPair
      (CaseTallWide $ \(Extent _ (EK.Square sh)) -> Left $ tall sh sh)
      (CaseTallWide Right)
      (CaseTallWide Left)
      (CaseTallWide $ \(Extent _ (EK.General h w)) ->
         if ge h w
            then Left $ tall h w
            else Right $ wide h w)


newtype GenSquare sh vert horiz =
   GenSquare {getGenSquare :: sh -> Extent vert horiz sh sh}

genSquare :: (C vert, C horiz) => sh -> Extent vert horiz sh sh
genSquare =
   getGenSquare $
   switchTagPair
      (GenSquare square)
      (GenSquare (\sh -> wide sh sh))
      (GenSquare (\sh -> tall sh sh))
      (GenSquare (\sh -> general sh sh))

newtype GenTall height width vert horiz =
   GenTall {
      getGenTall ::
         Extent vert Small height width -> Extent vert horiz height width
   }

generalizeTall :: (C vert, C horiz) =>
   Extent vert Small height width -> Extent vert horiz height width
generalizeTall =
   getGenTall $
   switchTagPair
      (GenTall id) (GenTall $ \(Extent _ (EK.Square s)) -> wide s s)
      (GenTall id) (GenTall $ \(Extent _ (EK.Tall h w)) -> general h w)

newtype GenWide height width vert horiz =
   GenWide {
      getGenWide ::
         Extent Small horiz height width -> Extent vert horiz height width
   }

generalizeWide :: (C vert, C horiz) =>
   Extent Small horiz height width -> Extent vert horiz height width
generalizeWide =
   getGenWide $
   switchTagPair
      (GenWide id)
      (GenWide id)
      (GenWide $ \(Extent _ (EK.Square s)) -> tall s s)
      (GenWide $ \(Extent _ (EK.Wide h w)) -> general h w)


newtype GenToTall height width vert horiz =
   GenToTall {
      getGenToTall ::
         Extent vert horiz height width -> Extent Big horiz height width
   }

genToTall :: (C vert, C horiz) =>
   Extent vert horiz height width -> Extent Big horiz height width
genToTall =
   getGenToTall $
   switchTagPair
      (GenToTall $ \(Extent _ (EK.Square s)) -> tall s s)
      (GenToTall $ \(Extent _ (EK.Wide h w)) -> general h w)
      (GenToTall id)
      (GenToTall id)


newtype GenToWide height width vert horiz =
   GenToWide {
      getGenToWide ::
         Extent vert horiz height width -> Extent vert Big height width
   }

genToWide :: (C vert, C horiz) =>
   Extent vert horiz height width -> Extent vert Big height width
genToWide =
   getGenToWide $
   switchTagPair
      (GenToWide $ \(Extent _ (EK.Square s)) -> wide s s)
      (GenToWide id)
      (GenToWide $ \(Extent _ (EK.Tall h w)) -> general h w)
      (GenToWide id)


squareSize :: Square sh -> sh
squareSize (Extent (Small,Small) (EK.Square sh)) = sh


newtype Accessor a height width vert horiz =
   Accessor {getAccessor :: Extent vert horiz height width -> a}

height :: (C vert, C horiz) => Extent vert horiz height width -> height
height =
   getAccessor $
   switchTagPair
      (Accessor (\(Extent _ (EK.Square s)) -> s))
      (Accessor (EK.wideHeight . extentDim))
      (Accessor (EK.tallHeight . extentDim))
      (Accessor (EK.generalHeight . extentDim))

width :: (C vert, C horiz) => Extent vert horiz height width -> width
width =
   getAccessor $
   switchTagPair
      (Accessor (\(Extent _ (EK.Square s)) -> s))
      (Accessor (EK.wideWidth . extentDim))
      (Accessor (EK.tallWidth . extentDim))
      (Accessor (EK.generalWidth . extentDim))


dimensions ::
   (C vert, C horiz) => Extent vert horiz height width -> (height,width)
dimensions x = (height x, width x)


toGeneral ::
   (C vert, C horiz) => Extent vert horiz height width -> General height width
toGeneral x = general (height x) (width x)

fromSquare :: (C vert, C horiz) => Square size -> Extent vert horiz size size
fromSquare = genSquare . squareSize

fromSquareLiberal :: (C vert, C horiz) =>
   Extent Small Small height width -> Extent vert horiz height width
fromSquareLiberal x@(Extent _ (EK.Square _)) = genSquare $ height x

squareFromGeneral ::
   (C vert, C horiz, Eq size) =>
   Extent vert horiz size size -> Square size
squareFromGeneral x =
   let size = height x
   in if size == width x
        then square size
        else error "Extent.squareFromGeneral: no square shape"


newtype Transpose height width vert horiz =
   Transpose {
      getTranspose ::
         Extent vert horiz height width ->
         Extent horiz vert width height
   }

transpose ::
   (C vert, C horiz) =>
   Extent vert horiz height width ->
   Extent horiz vert width height
transpose =
   getTranspose $
   switchTagPair
      (Transpose $ \(Extent o (EK.Square s)) -> Extent o (EK.Square s))
      (Transpose $ \(Extent o (EK.Wide h w)) -> Extent (swap o) (EK.Tall w h))
      (Transpose $ \(Extent o (EK.Tall h w)) -> Extent (swap o) (EK.Wide w h))
      (Transpose $ \(Extent o (EK.General h w)) -> Extent o (EK.General w h))


newtype Equal height width vert horiz =
   Equal {
      getEqual ::
         Extent vert horiz height width ->
         Extent vert horiz height width -> Bool
   }

instance
   (C vert, C horiz, Eq height, Eq width) =>
      Eq (Extent vert horiz height width) where
   (==) =
      getEqual $
      switchTagPair
         (Equal $ equating extentDim)
         (Equal $ equating extentDim)
         (Equal $ equating extentDim)
         (Equal $ equating extentDim)


instance
   (C vert, C horiz, Show height, Show width) =>
      Show (Extent vert horiz height width) where
   showsPrec prec =
      getAccessor $
      switchTagPair
         (Accessor $ showsPrecSquare prec)
         (Accessor $ showsPrecAny "Extent.wide" prec)
         (Accessor $ showsPrecAny "Extent.tall" prec)
         (Accessor $ showsPrecAny "Extent.general" prec)

showsPrecSquare ::
   (Show height) =>
   Int -> Extent Small Small height width -> ShowS
showsPrecSquare p x =
   showParen (p>10) $
   showString "Extent.square " . showsPrec 11 (height x)

showsPrecAny ::
   (C vert, C horiz, Show height, Show width) =>
   String -> Int -> Extent vert horiz height width -> ShowS
showsPrecAny name p x =
   showParen (p>10) $
   showString name .
   showString " " . showsPrec 11 (height x) .
   showString " " . showsPrec 11 (width x)


newtype Widen heightA widthA heightB widthB vert =
   Widen {
      getWiden ::
         Extent vert Big heightA widthA ->
         Extent vert Big heightB widthB
   }

widen ::
   (C vert) =>
   widthB -> Extent vert Big height widthA -> Extent vert Big height widthB
widen w =
   getWiden $
   switchTag
      (Widen (\(Extent o x) -> Extent o (x{EK.wideWidth = w})))
      (Widen (\(Extent o x) -> Extent o (x{EK.generalWidth = w})))

reduceWideHeight ::
   (C vert) =>
   heightB -> Extent vert Big heightA width -> Extent vert Big heightB width
reduceWideHeight h =
   getWiden $
   switchTag
      (Widen (\(Extent o x) -> Extent o (x{EK.wideHeight = h})))
      (Widen (\(Extent o x) -> Extent o (x{EK.generalHeight = h})))


newtype Adapt height width vert horiz =
   Adapt {
      getAdapt ::
         Extent vert horiz height width ->
         Extent vert horiz height width
   }

reduceConsistent ::
   (C vert, C horiz) =>
   height -> width ->
   Extent vert horiz height width -> Extent vert horiz height width
reduceConsistent h w =
   getAdapt $
   switchTagPair
      (Adapt $ \(Extent o (EK.Square _)) -> Extent o (EK.Square h))
      (Adapt $ \(Extent o (EK.Wide _ _)) -> Extent o (EK.Wide h w))
      (Adapt $ \(Extent o (EK.Tall _ _)) -> Extent o (EK.Tall h w))
      (Adapt $ \(Extent o (EK.General _ _)) -> Extent o (EK.General h w))


newtype Fuse height fuse width vert horiz =
   Fuse {
      getFuse ::
         Extent vert horiz height fuse ->
         Extent vert horiz fuse width ->
         Maybe (Extent vert horiz height width)
   }

fuse ::
   (C vert, C horiz, Eq fuse) =>
   Extent vert horiz height fuse ->
   Extent vert horiz fuse width ->
   Maybe (Extent vert horiz height width)
fuse =
   getFuse $
   switchTagPair
      (Fuse $
       \(Extent o (EK.Square s0)) (Extent _ (EK.Square s1)) ->
         toMaybe (s0==s1) $ Extent o (EK.Square s0))
      (Fuse $
       \(Extent o (EK.Wide h f0)) (Extent _ (EK.Wide f1 w)) ->
         toMaybe (f0==f1) $ Extent o (EK.Wide h w))
      (Fuse $
       \(Extent o (EK.Tall h f0)) (Extent _ (EK.Tall f1 w)) ->
         toMaybe (f0==f1) $ Extent o (EK.Tall h w))
      (Fuse $
       \(Extent o (EK.General h f0)) (Extent _ (EK.General f1 w)) ->
         toMaybe (f0==f1) $ Extent o (EK.General h w))


type family Multiply a b
type instance Multiply Small b = b
type instance Multiply Big   b = Big


data TagFact a = C a => TagFact

newtype MultiplyTagLaw b a =
   MultiplyTagLaw {
      getMultiplyTagLaw :: TagFact a -> TagFact b -> TagFact (Multiply a b)
   }

multiplyTagLaw :: TagFact a -> TagFact b -> TagFact (Multiply a b)
multiplyTagLaw a@TagFact =
   ($a) $ getMultiplyTagLaw $
   switchTag
      (MultiplyTagLaw $ flip const)
      (MultiplyTagLaw const)

heightFact :: (C vert) => Extent vert horiz height width -> TagFact vert
heightFact _ = TagFact

widthFact :: (C horiz) => Extent vert horiz height width -> TagFact horiz
widthFact _ = TagFact


newtype Unify height fuse width heightC widthC vertB horizB vertA horizA =
   Unify {
      getUnify ::
         Extent vertA horizA height fuse ->
         Extent vertB horizB fuse width ->
         Extent (Multiply vertA vertB) (Multiply horizA horizB) heightC widthC
   }

unifyLeft ::
   (C vertA, C horizA, C vertB, C horizB) =>
   Extent vertA horizA height fuse ->
   Extent vertB horizB fuse width ->
   Extent (Multiply vertA vertB) (Multiply horizA horizB) height fuse
unifyLeft =
   getUnify $
   switchTagPair
      (Unify $ const . fromSquareLiberal)
      (Unify $ const . generalizeWide)
      (Unify $ const . generalizeTall)
      (Unify $ const . toGeneral)

unifyRight ::
   (C vertA, C horizA, C vertB, C horizB) =>
   Extent vertA horizA height fuse ->
   Extent vertB horizB fuse width ->
   Extent (Multiply vertA vertB) (Multiply horizA horizB) fuse width
unifyRight =
   getUnify $
   switchTagPair
      (Unify $ const id)
      (Unify $ const genToWide)
      (Unify $ const genToTall)
      (Unify $ const toGeneral)


{-
Square  Square  -> Square
Square  Wide    -> Wide
Square  Tall    -> Tall
Square  General -> General
Wide    Square  -> Wide
Wide    Wide    -> Wide
Wide    Tall    -> General
Wide    General -> General
Tall    Square  -> Tall
Tall    Wide    -> General
Tall    Tall    -> Tall
Tall    General -> General
General Square  -> General
General Wide    -> General
General Tall    -> General
General General -> General

Small Small  Small Small -> Small Small
Small Small  Small Big   -> Small Big
Small Small  Big   Small -> Big   Small
Small Small  Big   Big   -> Big   Big
Small Big    Small Small -> Small Big
Small Big    Small Big   -> Small Big
Small Big    Big   Small -> Big   Big
Small Big    Big   Big   -> Big   Big
Big   Small  Small Small -> Big   Small
Big   Small  Small Big   -> Big   Big
Big   Small  Big   Small -> Big   Small
Big   Small  Big   Big   -> Big   Big
Big   Big    Small Small -> Big   Big
Big   Big    Small Big   -> Big   Big
Big   Big    Big   Small -> Big   Big
Big   Big    Big   Big   -> Big   Big
-}
