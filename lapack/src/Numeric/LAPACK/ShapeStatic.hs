{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.ShapeStatic where

import Numeric.LAPACK.Matrix.Shape.Private (UnaryProxy)

import qualified Data.FixedLength as FL

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)

import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num (integralFromProxy)
import Type.Base.Proxy (Proxy(Proxy))

import Foreign.Storable (Storable)

import Text.Printf (printf)


{- |
'ZeroBased' denotes a range starting at zero and has a certain length.
-}
newtype ZeroBased n = ZeroBased {zeroBasedSize :: UnaryProxy n}
   deriving (Eq, Show)

instance (Unary.Natural n) => Shape.C (ZeroBased n) where
   size = Shape.uncheckedSize
   uncheckedSize (ZeroBased len) = integralFromProxy len

instance (Unary.Natural n) => Shape.Indexed (ZeroBased n) where
   type Index (ZeroBased n) = FL.Index n
   indices _len = FL.toList FL.indices
   offset = Shape.uncheckedOffset
   uncheckedOffset _len = fromIntegral . FL.numFromIndex
   inBounds _len _ix = True

instance (Unary.Natural n) => Shape.InvIndexed (ZeroBased n) where
   -- could be implemented using new fixed-length-0.2.1:FL.indexFromNum
   indexFromOffset len k =
      case (0<=k, drop k $ Shape.indices len) of
         (True, i:_) -> i
         _ -> -- cf. comfort-array:Shape.errorIndexFromOffset
            error $
            printf "indexFromOffset (ShapeStatic.ZeroBased): index %d out of range" k


vector :: (Unary.Natural n, Storable a) => FL.T n a -> Array (ZeroBased n) a
vector = Array.fromList (ZeroBased Proxy) . FL.toList
