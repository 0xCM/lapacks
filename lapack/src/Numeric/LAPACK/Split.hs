module Numeric.LAPACK.Split where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Triangular.Private as TriPriv
import qualified Numeric.LAPACK.Matrix.Triangular.Basic as Tri
import qualified Numeric.LAPACK.Matrix.Private as Matrix
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Matrix.Triangular.Private (diagonalPointers)
import Numeric.LAPACK.Matrix.Triangular.Basic (UnitLower, Upper)
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor, ColumnMajor), transposeFromOrder,
          swapOnRowMajor, sideSwapFromOrder,
          Triangle, uploFromOrder, flipOrder)
import Numeric.LAPACK.Matrix.Private
         (Full, Transposition, transposeOrder,
          Conjugation(NonConjugated, Conjugated))
import Numeric.LAPACK.Linear.Private (solver, withInfo)
import Numeric.LAPACK.Scalar (zero, one)
import Numeric.LAPACK.Private (copyBlock, conjugateToTemp)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import System.IO.Unsafe (unsafePerformIO)

import Foreign.C.Types (CInt, CChar)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Ptr (Ptr)
import Foreign.Storable (poke)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)


type Split lower vert horiz height width =
      Array (MatrixShape.Split lower vert horiz height width)


determinantR ::
   (Extent.C vert, Shape.C height, Shape.C width, Class.Floating a) =>
   Split lower vert Extent.Small height width a -> a
determinantR (Array (MatrixShape.Split _ order extent) a) =
   let (height,width) = Extent.dimensions extent
       m = Shape.size height
       n = Shape.size width
       k = case order of RowMajor -> n; ColumnMajor -> m
   in unsafePerformIO $
      withForeignPtr a $ \aPtr ->
      Private.product (min m n) aPtr (k+1)

oddPermutation :: [CInt] -> Bool
oddPermutation = not . null . dropEven . filter id . zipWith (/=) [1..]

dropEven :: [a] -> [a]
dropEven (_:_:xs) = dropEven xs
dropEven xs = xs


extractTriangle ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Either lower Triangle ->
   Split lower vert horiz height width a ->
   Full vert horiz height width a
extractTriangle part (Array (MatrixShape.Split _ order extent) qr) =

   Array.unsafeCreate (MatrixShape.Full order extent) $ \rPtr -> do

   let (height,width) = Extent.dimensions extent
   let ((loup,m), (uplo,n)) =
         swapOnRowMajor order
            (('L', Shape.size height), ('U', Shape.size width))
   evalContT $ do
      loupPtr <- Call.char loup
      uploPtr <- Call.char uplo
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      qrPtr <- ContT $ withForeignPtr qr
      ldqrPtr <- Call.leadingDim m
      ldrPtr <- Call.leadingDim m
      zeroPtr <- Call.number zero
      onePtr <- Call.number one
      liftIO $
         case part of
            Left _ -> do
               LapackGen.lacpy loupPtr mPtr nPtr qrPtr ldqrPtr rPtr ldrPtr
               LapackGen.laset uploPtr mPtr nPtr zeroPtr onePtr rPtr ldrPtr
            Right _ -> do
               LapackGen.laset loupPtr mPtr nPtr zeroPtr zeroPtr rPtr ldrPtr
               LapackGen.lacpy uploPtr mPtr nPtr qrPtr ldqrPtr rPtr ldrPtr


wideExtractL ::
   (Extent.C horiz, Shape.C height, Shape.C width, Class.Floating a) =>
   Split lower Extent.Small horiz height width a -> UnitLower height a
wideExtractL =
   TriPriv.takeLower
      (MatrixShape.Unit,
       \order m lPtr -> mapM_ (flip poke one) $ diagonalPointers order m lPtr)
   .
   toFull

tallExtractR ::
   (Extent.C vert, Shape.C height, Shape.C width, Class.Floating a) =>
   Split lower vert Extent.Small height width a -> Upper width a
tallExtractR = Tri.takeUpper . toFull

toFull ::
   Split lower vert horiz height width a ->
   Full vert horiz height width a
toFull =
   Array.mapShape
      (\(MatrixShape.Split _ order extent) -> MatrixShape.Full order extent)


wideMultiplyL ::
   (Extent.C horizA, Extent.C vert, Extent.C horiz, Shape.C height, Eq height,
    Shape.C widthA, Shape.C widthB, Class.Floating a) =>
   Transposition ->
   Split Triangle Extent.Small horizA height widthA a ->
   Full vert horiz height widthB a ->
   Full vert horiz height widthB a
wideMultiplyL transposed a b =
   if MatrixShape.splitHeight (Array.shape a) == Matrix.height b
      then multiplyTriangular ('L','U') 'U' transposed a b
      else error "wideMultiplyL: height shapes mismatch"

tallMultiplyR ::
   (Extent.C vertA, Extent.C vert, Extent.C horiz, Shape.C height, Eq height,
    Shape.C heightA, Shape.C widthB, Class.Floating a) =>
   Transposition ->
   Split lower vertA Extent.Small heightA height a ->
   Full vert horiz height widthB a ->
   Full vert horiz height widthB a
tallMultiplyR transposed a b =
   if MatrixShape.splitWidth (Array.shape a) == Matrix.height b
      then multiplyTriangular ('U','L') 'N' transposed a b
      else error "wideMultiplyR: height shapes mismatch"

multiplyTriangular ::
   (Extent.C vertA, Extent.C horizA, Extent.C vertB, Extent.C horizB,
    Shape.C heightA, Shape.C widthA, Shape.C heightB, Shape.C widthB,
    Class.Floating a) =>
   (Char,Char) -> Char -> Transposition ->
   Split lower vertA horizA heightA widthA a ->
   Full vertB horizB heightB widthB a ->
   Full vertB horizB heightB widthB a
multiplyTriangular (normalPart,transposedPart) diag transposed
   (Array (MatrixShape.Split _ orderA extentA) a)
   (Array (MatrixShape.Full orderB extentB) b) =

   Array.unsafeCreate (MatrixShape.Full orderB extentB) $ \cPtr -> do

   let (heightA,widthA) = Extent.dimensions extentA
   let (heightB,widthB) = Extent.dimensions extentB
   let transOrderB = transposeOrder transposed orderB
   let ((uplo, transa), lda) =
         case orderA of
            RowMajor ->
               ((transposedPart, flipOrder transOrderB), Shape.size widthA)
            ColumnMajor ->
               ((normalPart, transOrderB), Shape.size heightA)
   let (side,(m,n)) =
         sideSwapFromOrder orderB (Shape.size heightB, Shape.size widthB)
   evalContT $ do
      sidePtr <- Call.char side
      uploPtr <- Call.char uplo
      transaPtr <- Call.char $ transposeFromOrder transa
      diagPtr <- Call.char diag
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim lda
      bPtr <- ContT $ withForeignPtr b
      ldcPtr <- Call.leadingDim m
      alphaPtr <- Call.number one
      liftIO $ do
         copyBlock (m*n) bPtr cPtr
         BlasGen.trmm sidePtr uploPtr transaPtr diagPtr
            mPtr nPtr alphaPtr aPtr ldaPtr cPtr ldcPtr


wideSolveL ::
   (Extent.C horizA, Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Shape.C nrhs, Class.Floating a) =>
   Transposition -> Conjugation ->
   Split Triangle Extent.Small horizA height width a ->
   Full vert horiz height nrhs a -> Full vert horiz height nrhs a
wideSolveL transposed conjugated
      (Array (MatrixShape.Split _ orderA extentA) a) =
   let heightA = Extent.height extentA
   in solver "Split.wideSolveL" heightA $ \n nPtr nrhsPtr xPtr ldxPtr -> do

      uploPtr <- Call.char $ uploFromOrder $ flipOrder orderA
      diagPtr <- Call.char 'U'
      let m = Shape.size heightA
      solveTriangular transposed conjugated orderA m n a
         uploPtr diagPtr nPtr nrhsPtr xPtr ldxPtr

tallSolveR ::
   (Extent.C vertA, Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Eq width, Shape.C nrhs, Class.Floating a) =>
   Transposition -> Conjugation ->
   Split lower vertA Extent.Small height width a ->
   Full vert horiz width nrhs a -> Full vert horiz width nrhs a
tallSolveR transposed conjugated
      (Array (MatrixShape.Split _ orderA extentA) a) =
   let (heightA,widthA) = Extent.dimensions extentA
   in solver "Split.tallSolveR" widthA $ \n nPtr nrhsPtr xPtr ldxPtr -> do

      uploPtr <- Call.char $ uploFromOrder orderA
      diagPtr <- Call.char 'N'
      let m = Shape.size heightA
      solveTriangular transposed conjugated orderA m n a
         uploPtr diagPtr nPtr nrhsPtr xPtr ldxPtr

solveTriangular ::
   Class.Floating a =>
   Transposition -> Conjugation ->
   Order -> Int -> Int -> ForeignPtr a ->
   Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt ->
   Ptr a -> Ptr CInt -> ContT r IO ()
solveTriangular transposed conjugated orderA m n a
   uploPtr diagPtr nPtr nrhsPtr xPtr ldxPtr = do
      let (trans, getA) =
            case (transposeOrder transposed orderA, conjugated) of
               (RowMajor, NonConjugated) -> ('T', ContT $ withForeignPtr a)
               (RowMajor, Conjugated) -> ('C', ContT $ withForeignPtr a)
               (ColumnMajor, NonConjugated) -> ('N', ContT $ withForeignPtr a)
               (ColumnMajor, Conjugated) -> ('N', conjugateToTemp (m*n) a)
      transPtr <- Call.char trans
      aPtr <- getA
      ldaPtr <- Call.leadingDim $ case orderA of ColumnMajor -> m; RowMajor -> n
      liftIO $
         withInfo "trtrs" $
            LapackGen.trtrs uploPtr transPtr diagPtr
               nPtr nrhsPtr aPtr ldaPtr xPtr ldxPtr
