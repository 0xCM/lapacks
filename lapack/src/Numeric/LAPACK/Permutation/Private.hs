module Numeric.LAPACK.Permutation.Private where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Split as Split
import Numeric.LAPACK.Matrix.Shape.Private (Order(RowMajor, ColumnMajor))
import Numeric.LAPACK.Matrix.Private
         (Full, Square, ZeroInt, Inversion(NonInverted, Inverted))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Format (Format(format))
import Numeric.LAPACK.Scalar (zero, one)
import Numeric.LAPACK.Private (fill, pointerSeq, copyBlock, copyToTemp)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import qualified Text.PrettyPrint.Boxes as TextBox

import qualified Foreign.Marshal.Array.Guarded as ForeignArray
import Foreign.Marshal.Array (advancePtr, copyArray)
import Foreign.C.Types (CInt)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Ptr (Ptr)
import Foreign.Storable (Storable, poke, peek, pokeElemOff, peekElemOff)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Monad (forM_)
import Control.Applicative ((<$>))

import Data.Bool.HT (if')


newtype Permutation sh = Permutation (Vector sh CInt)
   deriving (Show)

instance (Shape.C sh) => Format (Permutation sh) where
   format _fmt (Permutation perm) =
      let n = Shape.size $ Array.shape perm
      in TextBox.vcat TextBox.top $
         map (TextBox.hsep 1 TextBox.right . map TextBox.char) $
         map (\k -> (replicate (k-1) '.' ++ '1' : replicate (n-k) '.')) $
         map fromIntegral $ Array.toList perm


{-
We could use laswp if it would be available for CInt elements.
-}
{- |
The pivot array must be at most as long as @Shape.size sh@.
-}
fromPivots :: (Shape.C sh) =>
   Inversion -> sh -> Vector ZeroInt CInt -> Permutation sh
fromPivots inverted sh (Array (Shape.ZeroBased numIPiv) ipiv) =
   Permutation $
   if' (numIPiv > Shape.size sh)
      (error "Permutation.fromPivots: too many pivots") $
   Array.unsafeCreateWithSize sh $ \n permPtr ->
   withForeignPtr ipiv $ \ipivPtr -> do
      sequence_ $ take n $ zipWith poke (pointerSeq 1 permPtr) (iterate (1+) 1)
      let is =
            case inverted of
               Inverted -> tail $ iterate (subtract 1) numIPiv
               NonInverted -> iterate (1+) 0
      forM_ (take numIPiv is) $ \i ->
         swapElem permPtr i =<< peek1 ipivPtr i

swapElem :: (Storable a) => Ptr a -> Int -> Int -> IO ()
swapElem ptr i j = swap (advancePtr ptr i) (advancePtr ptr j)

swap :: (Storable a) => Ptr a -> Ptr a -> IO ()
swap ptr0 ptr1 = do
   a <- peek ptr0
   poke ptr0 =<< peek ptr1
   poke ptr1 a


toPivots :: (Shape.C sh) => Inversion -> Permutation sh -> Vector sh CInt
toPivots inverted (Permutation (Array sh perm)) =
   Array.unsafeCreateWithSize sh $ \n invPtr ->
   withForeignPtr perm $ \perm0Ptr ->
   ForeignArray.alloca n $ \permPtr -> do
      case inverted of
         Inverted -> do
            copyArray permPtr perm0Ptr n
            transposeIO n permPtr invPtr
         NonInverted -> do
            copyArray invPtr perm0Ptr n
            transposeIO n perm0Ptr permPtr
      forM_ (take n $ iterate (1+) 0) $ \i -> do
         j <- peek1 invPtr i
         k <- peek1 permPtr i
         poke1 permPtr j k
         poke1 invPtr k j


data Sign = Negative | Positive
   deriving (Eq, Show)

{-
We could also count the cycles of even number. This might be a little faster.
-}
determinant :: (Shape.C sh) => Permutation sh -> Sign
determinant =
   (\oddp -> if oddp then Negative else Positive) .
   Split.oddPermutation . Array.toList . toPivots NonInverted

numberFromSign :: (Class.Floating a) => Sign -> a
numberFromSign s =
   case s of
      Negative -> -1
      Positive -> 1


transpose :: (Shape.C sh) => Permutation sh -> Permutation sh
transpose (Permutation (Array shape perm)) =
   Permutation $
   Array.unsafeCreateWithSize shape $ \n dstPtr ->
   withForeignPtr perm $ \srcPtr ->
   transposeIO n srcPtr dstPtr

transposeIO :: Int -> Ptr CInt -> Ptr CInt -> IO ()
transposeIO n srcPtr dstPtr =
   forM_ (take n $ iterate (1+) 0) $ \i -> do
      j <- peek1 srcPtr i
      poke1 dstPtr j i


multiply :: (Shape.C sh, Eq sh) =>
   Permutation sh -> Permutation sh -> Permutation sh
multiply (Permutation (Array shape permA)) (Permutation (Array shapeB permB)) =
   if shape /= shapeB
      then error "Permutation.multiply: sizes mismatch"
      else
         Permutation $
         Array.unsafeCreateWithSize shape $ \n cPtr ->
         withForeignPtr permA $ \aPtr ->
         withForeignPtr permB $ \bPtr ->
         forM_ (take n $ iterate (1+) 0) $ \i ->
            poke1 cPtr i =<< peek1 bPtr =<< peek1 aPtr i


toMatrix :: (Shape.C sh, Class.Floating a) => Permutation sh -> Square sh a
toMatrix (Permutation (Array shape perm)) =
   Array.unsafeCreate (MatrixShape.square RowMajor shape) $ \aPtr ->
   withForeignPtr perm $ \permPtr -> do
      let n = Shape.size shape
      fill zero (n*n) aPtr
      forM_ (take n $ zip (iterate (1+) 0) (pointerSeq n aPtr)) $
         \(k,rowPtr) -> do
            i <- peek1 permPtr k
            pokeElemOff rowPtr i one


peek1 :: Ptr CInt -> Int -> IO Int
peek1 ptr i = subtract 1 . fromIntegral <$> peekElemOff ptr i

poke1 :: Ptr CInt -> Int -> Int -> IO ()
poke1 ptr i j = pokeElemOff ptr i (fromIntegral (j+1))


apply ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Class.Floating a) =>
   Bool -> Permutation height ->
   Full vert horiz height width a ->
   Full vert horiz height width a
apply inverted
      (Permutation (Array shapeP perm))
      (Array shape@(MatrixShape.Full order extent) a) =

   Array.unsafeCreateWithSize shape $ \blockSize bPtr -> do

   let (height,width) = Extent.dimensions extent
   Call.assert "Permutation.apply: heights mismatch" (height == shapeP)
   let m = Shape.size height
   let n = Shape.size width
   evalContT $ do
      fwdPtr <- Call.bool $ not inverted
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      kPtr <- copyToTemp n perm
      aPtr <- ContT $ withForeignPtr a
      liftIO $ do
         copyBlock blockSize aPtr bPtr
         case order of
            RowMajor -> LapackGen.lapmt fwdPtr nPtr mPtr bPtr mPtr kPtr
            ColumnMajor -> LapackGen.lapmr fwdPtr mPtr nPtr bPtr nPtr kPtr
