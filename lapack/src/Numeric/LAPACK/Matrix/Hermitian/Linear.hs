{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.Hermitian.Linear (
   solve,
   inverse,
   determinant,
   ) where

import Numeric.LAPACK.Matrix.Hermitian.Basic (Hermitian)
import Numeric.LAPACK.Matrix.Private (Full)

import qualified Numeric.LAPACK.Matrix.Symmetric.Private as Symmetric
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Hermitian.Private (Determinant(..))
import Numeric.LAPACK.Matrix.Private (Conjugation(Conjugated))
import Numeric.LAPACK.Scalar (RealOf, absoluteSquared)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import System.IO.Unsafe (unsafePerformIO)

import Foreign.Ptr (Ptr, castPtr)
import Foreign.Storable (peek)


solve ::
   (Extent.C vert, Extent.C horiz,
    Shape.C sh, Eq sh, Shape.C nrhs, Class.Floating a) =>
   Hermitian sh a ->
   Full vert horiz sh nrhs a -> Full vert horiz sh nrhs a
solve (Array (MatrixShape.Hermitian orderA shA) a) =
   Symmetric.solve "Hermitian.solve" Conjugated orderA shA a


inverse ::
   (Shape.C sh, Class.Floating a) => Hermitian sh a -> Hermitian sh a
inverse (Array shape@(MatrixShape.Hermitian order sh) a) =
   Array.unsafeCreateWithSize shape $
      Symmetric.inverse Conjugated order (Shape.size sh) a


determinant ::
   (Shape.C sh, Class.Floating a) => Hermitian sh a -> RealOf a
determinant =
   getDeterminant $
   Class.switchFloating
      (Determinant determinantAux) (Determinant determinantAux)
      (Determinant determinantAux) (Determinant determinantAux)

determinantAux ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Hermitian sh a -> ar
determinantAux (Array (MatrixShape.Hermitian order sh) a) =
   unsafePerformIO $
      Symmetric.determinant Conjugated
         peekBlockDeterminant order (Shape.size sh) a

peekBlockDeterminant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Ptr a, Maybe (Ptr a, Ptr a)) -> IO ar
peekBlockDeterminant (a0Ptr,ext) = do
   let peekReal = peek . castPtr
   a0 <- peekReal a0Ptr
   case ext of
      Nothing -> return a0
      Just (a1Ptr,bPtr) -> do
         a1 <- peekReal a1Ptr
         b <- peek bPtr
         return (a0*a1 - absoluteSquared b)
