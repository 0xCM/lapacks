{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.BandedHermitian.Eigen (
   values,
   decompose,
   ) where

import Numeric.LAPACK.Matrix.BandedHermitian.Basic (BandedHermitian)

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Private as Matrix
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

import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num (integralFromProxy)

import qualified Data.Array.Comfort.Storable.Unchecked.Monadic as ArrayIO
import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.C.Types (CInt, CChar)
import Foreign.Ptr (Ptr, nullPtr)
import Foreign.Storable (Storable)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)

import Data.Complex (Complex)


values ::
   (Unary.Natural offDiag, Shape.C sh, Class.Floating a) =>
   BandedHermitian offDiag sh a -> Vector sh (RealOf a)
values =
   runTakeDiagonal $
   Class.switchFloating
      (TakeDiagonal valuesAux) (TakeDiagonal valuesAux)
      (TakeDiagonal valuesAux) (TakeDiagonal valuesAux)

valuesAux ::
   (Unary.Natural offDiag, Shape.C sh,
    Class.Floating a, RealOf a ~ ar, Storable ar) =>
   BandedHermitian offDiag sh a -> Vector sh ar
valuesAux (Array (MatrixShape.BandedHermitian numOff order size) a) =
   Array.unsafeCreateWithSize size $ \n wPtr -> evalContT $ do
      let k = integralFromProxy numOff
      let lda = k+1
      jobzPtr <- Call.char 'N'
      uploPtr <- Call.char $ uploFromOrder order
      kPtr <- Call.cint k
      aPtr <- copyToTemp (n*lda) a
      ldaPtr <- Call.leadingDim lda
      let zPtr = nullPtr
      ldzPtr <- Call.leadingDim n
      liftIO $ withInfo eigenMsg "hbev" $
         hbev jobzPtr uploPtr n kPtr aPtr ldaPtr wPtr zPtr ldzPtr


decompose ::
   (Unary.Natural offDiag, Shape.C sh, Class.Floating a) =>
   BandedHermitian offDiag sh a -> (Matrix.Square sh a, Vector sh (RealOf a))
decompose =
   getDecompose $
   Class.switchFloating
      (Decompose decomposeAux) (Decompose decomposeAux)
      (Decompose decomposeAux) (Decompose decomposeAux)

type Decompose_ offDiag sh a =
      BandedHermitian offDiag sh a -> (Matrix.Square sh a, Vector sh (RealOf a))

newtype Decompose offDiag sh a =
   Decompose {getDecompose :: Decompose_ offDiag sh a}

decomposeAux ::
   (Unary.Natural offDiag, Shape.C sh,
    Class.Floating a, RealOf a ~ ar, Storable ar) =>
   Decompose_ offDiag sh a
decomposeAux (Array (MatrixShape.BandedHermitian numOff order size) a) =
   Array.unsafeCreateWithSizeAndResult (MatrixShape.square ColumnMajor size) $
      \_ zPtr ->
   ArrayIO.unsafeCreateWithSize size $ \n wPtr ->
   evalContT $ do
      let k = integralFromProxy numOff
      let lda = k+1
      jobzPtr <- Call.char 'V'
      uploPtr <- Call.char $ uploFromOrder order
      kPtr <- Call.cint k
      aPtr <- copyCondConjugateToTemp (order==RowMajor) (n*lda) a
      ldaPtr <- Call.leadingDim lda
      ldzPtr <- Call.leadingDim n
      liftIO $ withInfo eigenMsg "hbev" $
         hbev jobzPtr uploPtr n kPtr aPtr ldaPtr wPtr zPtr ldzPtr


type HBEV_ ar a =
   Ptr CChar -> Ptr CChar -> Int -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr ar ->
   Ptr a -> Ptr CInt -> Ptr CInt -> IO ()

newtype HBEV a = HBEV {getHBEV :: HBEV_ (RealOf a) a}

hbev :: Class.Floating a => HBEV_ (RealOf a) a
hbev =
   getHBEV $
   Class.switchFloating
      (HBEV sbevReal) (HBEV sbevReal) (HBEV hbevComplex) (HBEV hbevComplex)

sbevReal :: Class.Real a => HBEV_ a a
sbevReal jobzPtr uploPtr n kdPtr aPtr ldaPtr wPtr zPtr ldzPtr infoPtr =
   evalContT $ do
      nPtr <- Call.cint n
      workPtr <- Call.allocaArray (max 1 (3*n-2))
      liftIO $
         LapackReal.sbev jobzPtr uploPtr
            nPtr kdPtr aPtr ldaPtr wPtr zPtr ldzPtr workPtr infoPtr

hbevComplex :: Class.Real a => HBEV_ a (Complex a)
hbevComplex jobzPtr uploPtr n kdPtr aPtr ldaPtr wPtr zPtr ldzPtr infoPtr =
   evalContT $ do
      nPtr <- Call.cint n
      workPtr <- Call.allocaArray n
      rworkPtr <- Call.allocaArray (max 1 (3*n-2))
      liftIO $
         LapackComplex.hbev jobzPtr uploPtr
            nPtr kdPtr aPtr ldaPtr wPtr zPtr ldzPtr workPtr rworkPtr infoPtr
