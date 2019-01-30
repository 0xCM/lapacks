{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.Hermitian.Eigen (
   values,
   decompose,
   ) where

import Numeric.LAPACK.Matrix.Hermitian.Basic (Hermitian)
import Numeric.LAPACK.Matrix.Square.Basic (Square)

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import Numeric.LAPACK.Matrix.Hermitian.Private (TakeDiagonal(..))
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor,ColumnMajor), uploFromOrder)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf)
import Numeric.LAPACK.Private
         (copyToTemp, copyCondConjugateToTemp, withInfo, eigenMsg)

import qualified Numeric.LAPACK.FFI.Complex as LapackComplex
import qualified Numeric.LAPACK.FFI.Real as LapackReal
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked.Monadic as ArrayIO
import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))
import Data.Array.Comfort.Shape (triangleSize)

import Foreign.C.Types (CInt, CChar)
import Foreign.Ptr (Ptr, nullPtr)
import Foreign.Storable (Storable)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)

import Data.Complex (Complex)


values ::
   (Shape.C sh, Class.Floating a) =>
   Hermitian sh a -> Vector sh (RealOf a)
values =
   runTakeDiagonal $
   Class.switchFloating
      (TakeDiagonal valuesAux) (TakeDiagonal valuesAux)
      (TakeDiagonal valuesAux) (TakeDiagonal valuesAux)

valuesAux ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Storable ar) =>
   Hermitian sh a -> Vector sh ar
valuesAux (Array (MatrixShape.Hermitian order size) a) =
   Array.unsafeCreateWithSize size $ \n wPtr ->
   evalContT $ do
      jobzPtr <- Call.char 'N'
      uploPtr <- Call.char $ uploFromOrder order
      aPtr <- copyToTemp (triangleSize n) a
      let zPtr = nullPtr
      ldzPtr <- Call.leadingDim n
      liftIO $ withInfo eigenMsg "hpev" $
         hpev jobzPtr uploPtr n aPtr wPtr zPtr ldzPtr


decompose ::
   (Shape.C sh, Class.Floating a) =>
   Hermitian sh a -> (Square sh a, Vector sh (RealOf a))
decompose =
   getDecompose $
   Class.switchFloating
      (Decompose decomposeAux) (Decompose decomposeAux)
      (Decompose decomposeAux) (Decompose decomposeAux)

type Decompose_ sh a = Hermitian sh a -> (Square sh a, Vector sh (RealOf a))

newtype Decompose sh a = Decompose {getDecompose :: Decompose_ sh a}

decomposeAux ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Storable ar) =>
   Decompose_ sh a
decomposeAux (Array (MatrixShape.Hermitian order size) a) =
   Array.unsafeCreateWithSizeAndResult (MatrixShape.square ColumnMajor size) $
      \_ zPtr ->
   ArrayIO.unsafeCreateWithSize size $ \n wPtr ->
   evalContT $ do
      jobzPtr <- Call.char 'V'
      uploPtr <- Call.char $ uploFromOrder order
      aPtr <- copyCondConjugateToTemp (order==RowMajor) (triangleSize n) a
      ldzPtr <- Call.leadingDim n
      liftIO $ withInfo eigenMsg "hpev" $
         hpev jobzPtr uploPtr n aPtr wPtr zPtr ldzPtr


type HPEV_ ar a =
   Ptr CChar -> Ptr CChar -> Int -> Ptr a -> Ptr ar ->
   Ptr a -> Ptr CInt -> Ptr CInt -> IO ()

newtype HPEV a = HPEV {getHPEV :: HPEV_ (RealOf a) a}

hpev :: Class.Floating a => HPEV_ (RealOf a) a
hpev =
   getHPEV $
   Class.switchFloating
      (HPEV spevReal) (HPEV spevReal) (HPEV hpevComplex) (HPEV hpevComplex)

spevReal :: Class.Real a => HPEV_ a a
spevReal jobzPtr uploPtr n apPtr wPtr zPtr ldzPtr infoPtr =
   evalContT $ do
      nPtr <- Call.cint n
      workPtr <- Call.allocaArray (3*n)
      liftIO $
         LapackReal.spev
            jobzPtr uploPtr nPtr apPtr wPtr zPtr ldzPtr workPtr infoPtr

hpevComplex :: Class.Real a => HPEV_ a (Complex a)
hpevComplex jobzPtr uploPtr n apPtr wPtr zPtr ldzPtr infoPtr =
   evalContT $ do
      nPtr <- Call.cint n
      workPtr <- Call.allocaArray (max 1 (2*n-1))
      rworkPtr <- Call.allocaArray (max 1 (3*n-2))
      liftIO $
         LapackComplex.hpev
            jobzPtr uploPtr nPtr apPtr wPtr zPtr ldzPtr workPtr rworkPtr infoPtr
