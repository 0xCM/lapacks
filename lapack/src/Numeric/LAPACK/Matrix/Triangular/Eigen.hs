{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
module Numeric.LAPACK.Matrix.Triangular.Eigen (
   values,
   decompose,
   ) where

import qualified Numeric.LAPACK.Matrix.Triangular.Basic as Triangular
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import Numeric.LAPACK.Matrix.Triangular.Private
         (unpackZero, pack, unpackToTemp, fillTriangle,
          forPointers, rowMajorPointers)
import Numeric.LAPACK.Matrix.Triangular.Basic (Triangular)
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(ColumnMajor,RowMajor), caseLoUp, uploOrder, NonUnit(NonUnit))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (zero)
import Numeric.LAPACK.Private (lacgv, withInfo, errorCodeMsg)

import qualified Numeric.LAPACK.FFI.Complex as LapackComplex
import qualified Numeric.LAPACK.FFI.Real as LapackReal
import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked.Monadic as ArrayIO
import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))
import Data.Array.Comfort.Shape (triangleSize)

import Foreign.C.Types (CInt, CChar)
import Foreign.Ptr (Ptr, nullPtr)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)

import Data.Complex (Complex)
import Data.Tuple.HT (swap)


values ::
   (MatrixShape.DiagUpLo lo up, Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Vector sh a
values = Triangular.takeDiagonal


decompose ::
   (MatrixShape.DiagUpLo lo up, Shape.C sh, Class.Floating a) =>
   Triangular lo NonUnit up sh a ->
   (Triangular lo NonUnit up sh a, Vector sh a, Triangular lo NonUnit up sh a)
decompose a =
   let (vr,vl) =
         flip getDecompose a $
         MatrixShape.switchDiagUpLo
            (Decompose $
               (\eye -> (eye, Triangular.transpose eye)) .
               Triangular.relaxUnitDiagonal .
               Triangular.identity ColumnMajor .
               MatrixShape.triangularSize . Array.shape)
            (Decompose decomposeTriangular)
            (Decompose decomposeTriangular)
   in  (vr, values a, vl)

newtype Decompose sh a lo up =
   Decompose {
      getDecompose ::
         Triangular lo NonUnit up sh a ->
         (Triangular lo NonUnit up sh a, Triangular lo NonUnit up sh a)
   }

decomposeTriangular ::
   (MatrixShape.UpLo lo up, Shape.C sh, Class.Floating a) =>
   Triangular lo NonUnit up sh a ->
   (Triangular lo NonUnit up sh a, Triangular lo NonUnit up sh a)
decomposeTriangular (Array (MatrixShape.Triangular _diag uplo order sh) a) =
   let triShape ord =
         MatrixShape.Triangular NonUnit uplo (uploOrder uplo ord) sh
       n = Shape.size sh
       n2 = n*n
       triSize = triangleSize n

   in caseLoUp uplo id swap $
      Array.unsafeCreateWithSizeAndResult (triShape RowMajor) $ \_ vlpPtr ->
      ArrayIO.unsafeCreate (triShape ColumnMajor) $ \vrpPtr ->

   evalContT $ do
      sidePtr <- Call.char 'B'
      howManyPtr <- Call.char 'A'
      let selectPtr = nullPtr
      let unpk =
            case uploOrder uplo order of
               ColumnMajor -> unpackZero ColumnMajor
               RowMajor -> unpackZeroRowMajor
      aPtr <- unpackToTemp unpk n a
      ldaPtr <- Call.leadingDim n
      vlPtr <- Call.allocaArray n2
      vrPtr <- Call.allocaArray n2
      mmPtr <- Call.cint n
      mPtr <- Call.alloca
      liftIO $ withInfo errorCodeMsg "trevc" $
         trevc sidePtr howManyPtr selectPtr n
            aPtr ldaPtr vlPtr ldaPtr vrPtr ldaPtr mmPtr mPtr
      sizePtr <- Call.cint triSize
      incPtr <- Call.cint 1
      liftIO $ do
         pack ColumnMajor n vrPtr vrpPtr
         pack RowMajor n vlPtr vlpPtr
         lacgv sizePtr vlpPtr incPtr


unpackZeroRowMajor :: Class.Floating a => Int -> Ptr a -> Ptr a -> IO ()
unpackZeroRowMajor n packedPtr fullPtr = do
   fillTriangle zero RowMajor n fullPtr
   unpackRowMajor n packedPtr fullPtr

unpackRowMajor :: Class.Floating a => Int -> Ptr a -> Ptr a -> IO ()
unpackRowMajor n packedPtr fullPtr = evalContT $ do
   incxPtr <- Call.cint 1
   incyPtr <- Call.cint n
   liftIO $
      forPointers (rowMajorPointers n fullPtr packedPtr) $
            \nPtr (dstPtr,srcPtr) ->
         BlasGen.copy nPtr srcPtr incxPtr dstPtr incyPtr


type TREVC_ a =
   Ptr CChar -> Ptr CChar -> Ptr Bool ->
   Int -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt ->
   Ptr CInt -> Ptr CInt -> Ptr CInt -> IO ()

newtype TREVC a = TREVC {getTREVC :: TREVC_ a}

trevc :: Class.Floating a => TREVC_ a
trevc =
   getTREVC $
   Class.switchFloating
      (TREVC trevcReal) (TREVC trevcReal)
      (TREVC trevcComplex) (TREVC trevcComplex)

trevcReal :: Class.Real a => TREVC_ a
trevcReal sidePtr howmnyPtr selectPtr n
      tPtr ldtPtr vlPtr ldvlPtr vrPtr ldvrPtr mmPtr mPtr infoPtr =
   evalContT $ do
      nPtr <- Call.cint n
      workPtr <- Call.allocaArray (3*n)
      liftIO $
         LapackReal.trevc sidePtr howmnyPtr selectPtr nPtr
            tPtr ldtPtr vlPtr ldvlPtr vrPtr ldvrPtr mmPtr mPtr workPtr infoPtr

trevcComplex :: Class.Real a => TREVC_ (Complex a)
trevcComplex sidePtr howmnyPtr selectPtr n
      tPtr ldtPtr vlPtr ldvlPtr vrPtr ldvrPtr mmPtr mPtr infoPtr =
   evalContT $ do
      nPtr <- Call.cint n
      workPtr <- Call.allocaArray (2*n)
      rworkPtr <- Call.allocaArray n
      liftIO $
         LapackComplex.trevc sidePtr howmnyPtr selectPtr nPtr
            tPtr ldtPtr vlPtr ldvlPtr vrPtr ldvrPtr mmPtr mPtr
            workPtr rworkPtr infoPtr
