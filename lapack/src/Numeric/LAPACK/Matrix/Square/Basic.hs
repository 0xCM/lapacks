module Numeric.LAPACK.Matrix.Square.Basic (
   Square,
   size,
   toFull,
   toGeneral,
   fromGeneral,
   fromScalar,
   toScalar,
   fromList,
   autoFromList,

   transpose,
   adjoint,

   identity,
   identityFrom,
   identityFromWidth,
   identityFromHeight,
   diagonal,
   takeDiagonal,
   trace,

   multiply,
   square,
   power,
   ) where


import qualified Numeric.LAPACK.Matrix.Multiply as Mult
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as ExtentPriv
import qualified Numeric.LAPACK.Matrix.Extent as Extent
import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor, ColumnMajor), swapOnRowMajor)
import Numeric.LAPACK.Matrix.Private
         (Full, mapExtent,
          General, argGeneral, Square, argSquare, ZeroInt, zeroInt)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (zero, one)
import Numeric.LAPACK.Private (pokeCInt)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Storable (Storable, peek, poke)

import System.IO.Unsafe (unsafePerformIO)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)

import Data.Function.HT (powerAssociative)


size :: Square sh a -> sh
size = MatrixShape.fullHeight . Array.shape

toGeneral :: Square sh a -> General sh sh a
toGeneral = toFull

toFull ::
   (Extent.C vert, Extent.C horiz) => Square sh a -> Full vert horiz sh sh a
toFull = mapExtent Extent.fromSquare

fromGeneral :: (Eq sh) => General sh sh a -> Square sh a
fromGeneral = mapExtent (ExtentPriv.Map ExtentPriv.squareFromGeneral)


fromScalar :: (Storable a) => a -> Square () a
fromScalar a =
   Array.unsafeCreate (MatrixShape.square RowMajor ()) $ flip poke a

toScalar :: (Storable a) => Square () a -> a
toScalar = argSquare $ \_ () a ->
   unsafePerformIO $ withForeignPtr a peek

fromList :: (Shape.C sh, Storable a) => sh -> [a] -> Square sh a
fromList sh =
   Array.fromList (MatrixShape.square RowMajor sh)

autoFromList :: (Storable a) => [a] -> Square ZeroInt a
autoFromList xs =
   let n = length xs
       m = round $ sqrt (fromIntegral n :: Double)
   in if n == m*m
        then fromList (zeroInt m) xs
        else error "Square.autoFromList: no quadratic number of elements"


transpose :: Square sh a -> Square sh a
transpose = Array.mapShape MatrixShape.transpose

{- |
conjugate transpose
-}
adjoint :: (Shape.C sh, Class.Floating a) => Square sh a -> Square sh a
adjoint = transpose . Vector.conjugate


identity :: (Shape.C sh, Class.Floating a) => sh -> Square sh a
identity = identityOrder ColumnMajor

identityFrom :: (Shape.C sh, Class.Floating a) => Square sh a -> Square sh a
identityFrom = argSquare $ \order sh _ -> identityOrder order sh

identityFromWidth ::
   (Shape.C height, Shape.C width, Class.Floating a) =>
   General height width a -> Square width a
identityFromWidth =
   argGeneral $ \order _ width _ -> identityOrder order width

identityFromHeight ::
   (Shape.C height, Shape.C width, Class.Floating a) =>
   General height width a -> Square height a
identityFromHeight =
   argGeneral $ \order height _ _ -> identityOrder order height

identityOrder, _identityOrder ::
   (Shape.C sh, Class.Floating a) => Order -> sh -> Square sh a
identityOrder order sh =
   Array.unsafeCreate (MatrixShape.square order sh) $ \aPtr ->
   evalContT $ do
      uploPtr <- Call.char 'A'
      nPtr <- Call.cint $ Shape.size sh
      alphaPtr <- Call.number zero
      betaPtr <- Call.number one
      liftIO $ LapackGen.laset uploPtr nPtr nPtr alphaPtr betaPtr aPtr nPtr

_identityOrder order sh =
   Array.unsafeCreateWithSize (MatrixShape.square order sh) $ \blockSize yPtr ->
   evalContT $ do
      nPtr <- Call.alloca
      xPtr <- Call.number zero
      incxPtr <- Call.cint 0
      incyPtr <- Call.cint 1
      liftIO $ do
         pokeCInt nPtr blockSize
         BlasGen.copy nPtr xPtr incxPtr yPtr incyPtr
         let n = fromIntegral $ Shape.size sh
         poke nPtr n
         poke xPtr one
         poke incyPtr (n+1)
         BlasGen.copy nPtr xPtr incxPtr yPtr incyPtr

diagonal :: (Shape.C sh, Class.Floating a) => Vector sh a -> Square sh a
diagonal (Array sh x) =
   Array.unsafeCreateWithSize (MatrixShape.square ColumnMajor sh) $
      \blockSize yPtr ->
   evalContT $ do
      nPtr <- Call.alloca
      xPtr <- ContT $ withForeignPtr x
      zPtr <- Call.number zero
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 1
      inczPtr <- Call.cint 0
      liftIO $ do
         pokeCInt nPtr blockSize
         BlasGen.copy nPtr zPtr inczPtr yPtr incyPtr
         let n = fromIntegral $ Shape.size sh
         poke nPtr n
         poke incyPtr (n+1)
         BlasGen.copy nPtr xPtr incxPtr yPtr incyPtr

takeDiagonal :: (Shape.C sh, Class.Floating a) => Square sh a -> Vector sh a
takeDiagonal = argSquare $ \_ sh x ->
   Array.unsafeCreateWithSize sh $ \n yPtr -> evalContT $ do
      nPtr <- Call.cint n
      xPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint (n+1)
      incyPtr <- Call.cint 1
      liftIO $ BlasGen.copy nPtr xPtr incxPtr yPtr incyPtr

trace :: (Shape.C sh, Class.Floating a) => Square sh a -> a
trace = argSquare $ \_ sh x -> unsafePerformIO $ do
   let n = Shape.size sh
   withForeignPtr x $ \xPtr -> Private.sum n xPtr (n+1)


multiply ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Square sh a -> Square sh a -> Square sh a
multiply = Mult.multiply

square :: (Shape.C sh, Class.Floating a) => Square sh a -> Square sh a
square a = multiplyCommutativeUnchecked a a

power ::
   (Shape.C sh, Class.Floating a) =>
   Integer -> Square sh a -> Square sh a
power n a =
   powerAssociative multiplyCommutativeUnchecked (identityFrom a) a n

{-
orderA and orderB must be equal but this is not checked.
-}
multiplyCommutativeUnchecked ::
   (Shape.C sh, Class.Floating a) =>
   Square sh a -> Square sh a -> Square sh a
multiplyCommutativeUnchecked
   (Array shape@(MatrixShape.Full order extent) a)
   (Array _ b) =
      Array.unsafeCreate shape $ \cPtr ->
   let n = Shape.size $ Extent.height extent
       (at,bt) = swapOnRowMajor order (a,b)
   in  Private.multiplyMatrix ColumnMajor ColumnMajor n n n at bt cPtr
