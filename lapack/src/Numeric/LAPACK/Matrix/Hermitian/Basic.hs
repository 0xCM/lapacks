{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.Hermitian.Basic (
   Hermitian,
   Transposition(..),
   fromList,
   autoFromList,
   identity,
   diagonal,
   takeDiagonal,

   multiplyVector,
   square,
   multiplyFull,
   outer,
   sumRank1, sumRank1NonEmpty,
   sumRank2, sumRank2NonEmpty,

   toSquare,
   covariance,
   addAdjoint,
   ) where

import qualified Numeric.LAPACK.Matrix.Symmetric.Private as Symmetric
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Hermitian.Private (Diagonal(..), TakeDiagonal(..))
import Numeric.LAPACK.Matrix.Triangular.Private
         (forPointers, pack, unpack, unpackToTemp,
          diagonalPointers, diagonalPointerPairs,
          rowMajorPointers, columnMajorPointers)
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor,ColumnMajor), flipOrder, sideSwapFromOrder,
          uploFromOrder)
import Numeric.LAPACK.Matrix.Private
         (Full, General, argGeneral, Square, argSquare, ZeroInt, zeroInt,
          Transposition(NonTransposed, Transposed), transposeOrder,
          Conjugation(Conjugated))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, zero, one, fromReal, realPart)
import Numeric.LAPACK.Private
         (fill, lacgv, copyConjugate, conjugateToTemp, condConjugateToTemp)

import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.BLAS.FFI.Complex as BlasComplex
import qualified Numeric.BLAS.FFI.Real as BlasReal
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.C.Types (CInt, CChar)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Ptr (Ptr)
import Foreign.Storable (Storable, poke, peek)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Monad (when)

import qualified Data.NonEmpty as NonEmpty
import Data.Foldable (forM_)


type Hermitian sh = Array (MatrixShape.Hermitian sh)


fromList :: (Shape.C sh, Storable a) => Order -> sh -> [a] -> Hermitian sh a
fromList order sh =
   Array.fromList (MatrixShape.Hermitian order sh)

autoFromList :: (Storable a) => Order -> [a] -> Hermitian ZeroInt a
autoFromList order xs =
   fromList order
      (zeroInt $ MatrixShape.triangleExtent "Hermitian.autoFromList" $
       length xs)
      xs


identity :: (Shape.C sh, Class.Floating a) => Order -> sh -> Hermitian sh a
identity order sh =
   Array.unsafeCreateWithSize (MatrixShape.Hermitian order sh) $
      \triSize aPtr -> do
   fill zero triSize aPtr
   mapM_ (flip poke one) $ diagonalPointers order (Shape.size sh) aPtr

diagonal ::
   (Shape.C sh, Class.Floating a) =>
   Order -> Vector sh (RealOf a) -> Hermitian sh a
diagonal order =
   runDiagonal $
   Class.switchFloating
      (Diagonal $ diagonalAux order) (Diagonal $ diagonalAux order)
      (Diagonal $ diagonalAux order) (Diagonal $ diagonalAux order)

diagonalAux ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Storable ar) =>
   Order -> Vector sh ar -> Hermitian sh a
diagonalAux order (Array sh x) =
   Array.unsafeCreateWithSize (MatrixShape.Hermitian order sh) $
      \triSize aPtr -> do
   fill zero triSize aPtr
   withForeignPtr x $ \xPtr ->
      forM_ (diagonalPointerPairs order (Shape.size sh) xPtr aPtr) $
         \(srcPtr,dstPtr) -> poke dstPtr . fromReal =<< peek srcPtr


takeDiagonal ::
   (Shape.C sh, Class.Floating a) =>
   Hermitian sh a -> Vector sh (RealOf a)
takeDiagonal =
   runTakeDiagonal $
   Class.switchFloating
      (TakeDiagonal takeDiagonalAux) (TakeDiagonal takeDiagonalAux)
      (TakeDiagonal takeDiagonalAux) (TakeDiagonal takeDiagonalAux)

takeDiagonalAux ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Storable ar) =>
   Hermitian sh a -> Vector sh ar
takeDiagonalAux (Array (MatrixShape.Hermitian order sh) a) =
   Array.unsafeCreateWithSize sh $ \n xPtr ->
   withForeignPtr a $ \aPtr ->
      forM_ (diagonalPointerPairs order n xPtr aPtr) $
         \(dstPtr,srcPtr) -> poke dstPtr . realPart =<< peek srcPtr


multiplyVector ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Transposition -> Hermitian sh a -> Vector sh a -> Vector sh a
multiplyVector transposed
   (Array (MatrixShape.Hermitian order shA) a) (Array shX x) =
      Array.unsafeCreateWithSize shX $ \n yPtr -> do
   Call.assert "Hermitian.multiplyVector: width shapes mismatch" (shA == shX)
   evalContT $ do
      let conj = transposeOrder transposed order == RowMajor
      uploPtr <- Call.char $ uploFromOrder order
      nPtr <- Call.cint n
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      xPtr <- condConjugateToTemp conj n x
      incxPtr <- Call.cint 1
      betaPtr <- Call.number zero
      incyPtr <- Call.cint 1
      liftIO $ do
         BlasGen.hpmv
            uploPtr nPtr alphaPtr aPtr xPtr incxPtr betaPtr yPtr incyPtr
         when conj $ lacgv nPtr yPtr incyPtr


square ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Hermitian sh a -> Hermitian sh a
square (Array shape@(MatrixShape.Hermitian order sh) a) =
   Array.unsafeCreate shape $
      Symmetric.square Conjugated order (Shape.size sh) a


multiplyFull ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width,
    Class.Floating a) =>
   Transposition -> Hermitian height a ->
   Full vert horiz height width a ->
   Full vert horiz height width a
multiplyFull transposed
   (Array        (MatrixShape.Hermitian orderA shA) a)
   (Array shapeB@(MatrixShape.Full orderB extentB) b) =
      Array.unsafeCreate shapeB $ \cPtr -> do
   let (height,width) = Extent.dimensions extentB
   Call.assert "Hermitian.multiplyFull: shapes mismatch" (shA == height)
   let m0 = Shape.size height
   let n0 = Shape.size width
   let size = m0*m0
   evalContT $ do
      let (side,(m,n)) = sideSwapFromOrder orderB (m0,n0)
      sidePtr <- Call.char side
      uploPtr <- Call.char $ uploFromOrder orderA
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.number one
      aPtr <- unpackToTemp (unpack orderA) m0 a
      ldaPtr <- Call.leadingDim m0
      incaPtr <- Call.cint 1
      sizePtr <- Call.cint size
      bPtr <- ContT $ withForeignPtr b
      ldbPtr <- Call.leadingDim m
      betaPtr <- Call.number zero
      ldcPtr <- Call.leadingDim m
      liftIO $ do
         when (transposeOrder transposed orderA /= orderB) $
            lacgv sizePtr aPtr incaPtr
         BlasGen.hemm sidePtr uploPtr
            mPtr nPtr alphaPtr aPtr ldaPtr
            bPtr ldbPtr betaPtr cPtr ldcPtr



withConjBuffer ::
   (Shape.C sh, Class.Floating a) =>
   Order -> sh -> Int -> Ptr a ->
   (Ptr CChar -> Ptr CInt -> Ptr CInt -> IO ()) -> ContT r IO ()
withConjBuffer order sh triSize aPtr act = do
   uploPtr <- Call.char $ uploFromOrder order
   nPtr <- Call.cint $ Shape.size sh
   incxPtr <- Call.cint 1
   sizePtr <- Call.cint triSize
   liftIO $ do
      fill zero triSize aPtr
      act uploPtr nPtr incxPtr
      case order of
         RowMajor -> lacgv sizePtr aPtr incxPtr
         ColumnMajor -> return ()

outer ::
   (Shape.C sh, Class.Floating a) => Order -> Vector sh a -> Hermitian sh a
outer order =
   getMap $
   Class.switchFloating
      (Map $ outerAux order) (Map $ outerAux order)
      (Map $ outerAux order) (Map $ outerAux order)

outerAux ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Order -> Vector sh a -> Hermitian sh a
outerAux order (Array sh x) =
   Array.unsafeCreateWithSize (MatrixShape.Hermitian order sh) $
      \triSize aPtr ->
   evalContT $ do
      alphaPtr <- Call.real one
      xPtr <- ContT $ withForeignPtr x
      withConjBuffer order sh triSize aPtr $ \uploPtr nPtr incxPtr ->
         hpr uploPtr nPtr alphaPtr xPtr incxPtr aPtr


sumRank1 ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Order -> sh -> [(RealOf a, Vector sh a)] -> Hermitian sh a
sumRank1 =
   getSumRank1 $
   Class.switchFloating
      (SumRank1 sumRank1Aux) (SumRank1 sumRank1Aux)
      (SumRank1 sumRank1Aux) (SumRank1 sumRank1Aux)

type SumRank1_ sh ar a = Order -> sh -> [(ar, Vector sh a)] -> Hermitian sh a

newtype SumRank1 sh a = SumRank1 {getSumRank1 :: SumRank1_ sh (RealOf a) a}

sumRank1Aux ::
   (Shape.C sh, Eq sh, Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   SumRank1_ sh ar a
sumRank1Aux order sh xs =
   Array.unsafeCreateWithSize (MatrixShape.Hermitian order sh) $
      \triSize aPtr ->
   evalContT $ do
      alphaPtr <- Call.alloca
      withConjBuffer order sh triSize aPtr $ \uploPtr nPtr incxPtr ->
         forM_ xs $ \(alpha, Array shX x) ->
            withForeignPtr x $ \xPtr -> do
               Call.assert
                  "Hermitian.sumRank1: non-matching vector size" (sh==shX)
               poke alphaPtr alpha
               hpr uploPtr nPtr alphaPtr xPtr incxPtr aPtr


sumRank1NonEmpty ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Order -> NonEmpty.T [] (RealOf a, Vector sh a) -> Hermitian sh a
sumRank1NonEmpty order (NonEmpty.Cons x xs) =
   sumRank1 order (Array.shape $ snd x) (x:xs)


type HPR_ a =
   Ptr CChar -> Ptr CInt ->
   Ptr (RealOf a) -> Ptr a -> Ptr CInt -> Ptr a -> IO ()

newtype HPR a = HPR {getHPR :: HPR_ a}

hpr :: Class.Floating a => HPR_ a
hpr =
   getHPR $
   Class.switchFloating
      (HPR BlasReal.spr) (HPR BlasReal.spr)
      (HPR BlasComplex.hpr) (HPR BlasComplex.hpr)


sumRank2 ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Order -> sh -> [(a, (Vector sh a, Vector sh a))] -> Hermitian sh a
sumRank2 order sh xys =
   Array.unsafeCreateWithSize (MatrixShape.Hermitian order sh) $
      \triSize aPtr ->
   evalContT $ do
      alphaPtr <- Call.alloca
      withConjBuffer order sh triSize aPtr $ \uploPtr nPtr incPtr ->
         forM_ xys $ \(alpha, (Array shX x, Array shY y)) ->
            withForeignPtr x $ \xPtr ->
            withForeignPtr y $ \yPtr -> do
               Call.assert
                  "Hermitian.sumRank2: non-matching x vector size" (sh==shX)
               Call.assert
                  "Hermitian.sumRank2: non-matching y vector size" (sh==shY)
               poke alphaPtr alpha
               BlasGen.hpr2 uploPtr nPtr alphaPtr xPtr incPtr yPtr incPtr aPtr

sumRank2NonEmpty ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Order -> NonEmpty.T [] (a, (Vector sh a, Vector sh a)) -> Hermitian sh a
sumRank2NonEmpty order (NonEmpty.Cons xy xys) =
   sumRank2 order (Array.shape $ fst $ snd xy) (xy:xys)


{-
It is not strictly necessary to keep the 'order'.
It would be neither more complicated nor less efficient
to change the order via the conversion.
-}
toSquare, _toSquare ::
   (Shape.C sh, Class.Floating a) => Hermitian sh a -> Square sh a
_toSquare (Array (MatrixShape.Hermitian order sh) a) =
      Array.unsafeCreate (MatrixShape.square order sh) $ \bPtr ->
   evalContT $ do
      let n = Shape.size sh
      aPtr <- ContT $ withForeignPtr a
      conjPtr <- conjugateToTemp (Shape.triangleSize n) a
      liftIO $ do
         unpack (flipOrder order) n conjPtr bPtr -- wrong
         unpack order n aPtr bPtr

toSquare (Array (MatrixShape.Hermitian order sh) a) =
      Array.unsafeCreate (MatrixShape.square order sh) $ \bPtr ->
   withForeignPtr a $ \aPtr ->
      Symmetric.unpack Conjugated order (Shape.size sh) aPtr bPtr


{- |
A^H * A
-}
covariance ::
   (Shape.C height, Shape.C width, Eq width, Class.Floating a) =>
   General height width a -> Hermitian width a
covariance =
   getMap $
   Class.switchFloating
      (Map covarianceAux) (Map covarianceAux)
      (Map covarianceAux) (Map covarianceAux)

newtype Map f g a = Map {getMap :: f a -> g a}

covarianceAux ::
   (Shape.C height, Shape.C width, Eq width,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General height width a -> Hermitian width a
covarianceAux = argGeneral $ \order height width a ->
   Array.unsafeCreate (MatrixShape.Hermitian order width) $ \bPtr -> do

   let n = Shape.size width
   let k = Shape.size height
   evalContT $ do
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      betaPtr <- Call.number zero
      cPtr <- Call.allocaArray (n*n)
      ldcPtr <- Call.leadingDim n

      case order of
         ColumnMajor -> do
            uploPtr <- Call.char 'U'
            transPtr <- Call.char 'C'
            ldaPtr <- Call.leadingDim k
            liftIO $ do
               herk uploPtr transPtr
                  nPtr kPtr alphaPtr aPtr ldaPtr betaPtr cPtr ldcPtr
               pack ColumnMajor n cPtr bPtr

         RowMajor -> do
            uploPtr <- Call.char 'L'
            transPtr <- Call.char 'N'
            ldaPtr <- Call.leadingDim n
            liftIO $ do
               herk uploPtr transPtr
                  nPtr kPtr alphaPtr aPtr ldaPtr betaPtr cPtr ldcPtr
               pack RowMajor n cPtr bPtr


type HERK_ a =
   Ptr CChar -> Ptr CChar -> Ptr CInt -> Ptr CInt -> Ptr (RealOf a) -> Ptr a ->
   Ptr CInt -> Ptr (RealOf a) -> Ptr a -> Ptr CInt -> IO ()

newtype HERK a = HERK {getHERK :: HERK_ a}

herk :: Class.Floating a => HERK_ a
herk =
   getHERK $
   Class.switchFloating
      (HERK BlasReal.syrk)
      (HERK BlasReal.syrk)
      (HERK BlasComplex.herk)
      (HERK BlasComplex.herk)


{- |
A^H + A
-}
addAdjoint, _addAdjoint ::
   (Shape.C sh, Class.Floating a) => Square sh a -> Hermitian sh a
_addAdjoint =
   argSquare $ \order sh a ->
      Array.unsafeCreateWithSize (MatrixShape.Hermitian order sh) $ \bSize bPtr -> do
   let n = Shape.size sh
   evalContT $ do
      alphaPtr <- Call.number one
      incxPtr <- Call.cint 1
      aPtr <- ContT $ withForeignPtr a
      sizePtr <- Call.cint bSize
      conjPtr <- Call.allocaArray bSize
      liftIO $ do
         pack order n aPtr bPtr
         pack (flipOrder order) n aPtr conjPtr -- wrong
         lacgv sizePtr conjPtr incxPtr
         BlasGen.axpy sizePtr alphaPtr conjPtr incxPtr bPtr incxPtr

addAdjoint =
   argSquare $ \order sh a ->
      Array.unsafeCreate (MatrixShape.Hermitian order sh) $ \bPtr -> do
   let n = Shape.size sh
   evalContT $ do
      alphaPtr <- Call.number one
      incxPtr <- Call.cint 1
      incnPtr <- Call.cint n
      aPtr <- ContT $ withForeignPtr a
      liftIO $ case order of
         RowMajor ->
            forPointers (rowMajorPointers n aPtr bPtr) $
               \nPtr (srcPtr,dstPtr) -> do
            copyConjugate nPtr srcPtr incnPtr dstPtr incxPtr
            BlasGen.axpy nPtr alphaPtr srcPtr incxPtr dstPtr incxPtr
         ColumnMajor ->
            forPointers (columnMajorPointers n aPtr bPtr) $
               \nPtr ((srcRowPtr,srcColumnPtr),dstPtr) -> do
            copyConjugate nPtr srcRowPtr incnPtr dstPtr incxPtr
            BlasGen.axpy nPtr alphaPtr srcColumnPtr incxPtr dstPtr incxPtr


_pack :: Class.Floating a => Order -> Int -> Ptr a -> Ptr a -> IO ()
_pack order n fullPtr packedPtr =
   evalContT $ do
      incxPtr <- Call.cint 1
      liftIO $
         case order of
            ColumnMajor ->
               forPointers (columnMajorPointers n fullPtr packedPtr) $
                  \nPtr ((_,srcPtr),dstPtr) ->
                     BlasGen.copy nPtr srcPtr incxPtr dstPtr incxPtr
            RowMajor ->
               forPointers (rowMajorPointers n fullPtr packedPtr) $
                  \nPtr (srcPtr,dstPtr) ->
                     BlasGen.copy nPtr srcPtr incxPtr dstPtr incxPtr
