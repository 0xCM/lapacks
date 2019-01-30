{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE GADTs #-}
module Numeric.LAPACK.Matrix.BandedHermitian.Basic (
   BandedHermitian,
   Transposition(..),
   fromList,
   identity,
   diagonal,
   takeDiagonal,
   toHermitian,
   toBanded,
   multiplyVector,
   multiplyFull,
   covariance,
   sumRank1,
   ) where

import qualified Numeric.LAPACK.ShapeStatic as ShapeStatic
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent as Extent
import qualified Numeric.LAPACK.Matrix.Banded.Basic as Banded
import qualified Numeric.LAPACK.Matrix.Triangular.Private as TriangularPriv
import qualified Numeric.LAPACK.Matrix.Private as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix.Hermitian.Private (TakeDiagonal(..))
import Numeric.LAPACK.Matrix.Hermitian.Basic (Hermitian)
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor,ColumnMajor), flipOrder, uploFromOrder,
          UnaryProxy, natFromProxy)
import Numeric.LAPACK.Matrix.Private
         (Transposition(NonTransposed, Transposed), transposeOrder)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, zero, one)
import Numeric.LAPACK.Private
         (fill, lacgv, copyConjugate, condConjugateToTemp,
          pointerSeq, pokeCInt, copySubMatrix)

import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.BLAS.FFI.Complex as BlasComplex
import qualified Numeric.BLAS.FFI.Real as BlasReal
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary.Literal as TypeNum
import qualified Type.Data.Num.Unary.Proof as Proof
import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num.Unary ((:+:))
import Type.Data.Num (integralFromProxy)
import Type.Base.Proxy (Proxy(Proxy))

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.Marshal.Array (advancePtr)
import Foreign.C.Types (CInt, CChar)
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Ptr (Ptr, castPtr)
import Foreign.Storable (Storable, poke, peek, peekElemOff)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Monad (when)

import Data.Foldable (for_)
import Data.Tuple.HT (mapPair)

import Data.Complex (Complex, conjugate)


type BandedHermitian offDiag size =
      Array (MatrixShape.BandedHermitian offDiag size)

type Diagonal size = BandedHermitian TypeNum.U0 size


fromList ::
   (Unary.Natural offDiag, Shape.C size, Storable a) =>
   UnaryProxy offDiag -> Order -> size -> [a] ->
   BandedHermitian offDiag size a
fromList numOff order size =
   Array.fromList (MatrixShape.BandedHermitian numOff order size)

identity ::
   (Shape.C sh, Class.Floating a) => sh -> Diagonal sh a
identity sh =
   Array.mapShape (MatrixShape.BandedHermitian Proxy ColumnMajor) $
   Vector.constant sh one

diagonal ::
   (Shape.C sh, Class.Floating a) => Vector sh (RealOf a) -> Diagonal sh a
diagonal =
   Array.mapShape (MatrixShape.BandedHermitian Proxy ColumnMajor) .
   Vector.fromReal

takeDiagonal ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   BandedHermitian offDiag size a -> Vector size (RealOf a)
takeDiagonal =
   runTakeDiagonal $
   Class.switchFloating
      (TakeDiagonal $ takeDiagonalAux 1) (TakeDiagonal $ takeDiagonalAux 1)
      (TakeDiagonal $ takeDiagonalAux 2) (TakeDiagonal $ takeDiagonalAux 2)

takeDiagonalAux ::
   (Unary.Natural offDiag, Shape.C size,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Int -> BandedHermitian offDiag size a -> Vector size ar
takeDiagonalAux dim (Array (MatrixShape.BandedHermitian numOff order size) a) =
   let k = integralFromProxy numOff
   in Array.unsafeCreateWithSize size $ \n yPtr -> evalContT $ do
         nPtr <- Call.cint n
         aPtr <- ContT $ withForeignPtr a
         let xPtr =
               castPtr $ advancePtr aPtr $
               case order of
                  RowMajor -> 0
                  ColumnMajor -> k
         incxPtr <- Call.cint (dim * (k+1))
         incyPtr <- Call.cint 1
         liftIO $ BlasGen.copy nPtr xPtr incxPtr yPtr incyPtr


toHermitian ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   BandedHermitian offDiag size a -> Hermitian size a
toHermitian (Array (MatrixShape.BandedHermitian numOff order size) a) =
   Array.unsafeCreateWithSize (MatrixShape.Hermitian order size) $
   TriangularPriv.fromBanded
      (integralFromProxy numOff) order (Shape.size size) a


toBanded ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   BandedHermitian offDiag size a ->
   Banded.Square offDiag offDiag size a
toBanded (Array (MatrixShape.BandedHermitian numOff order sh) a) =
   Array.unsafeCreate
      (MatrixShape.Banded (numOff,numOff) order (Extent.square sh)) $ \bPtr ->
   withForeignPtr a $ \aPtr ->
      case order of
         ColumnMajor -> toBandedColumnMajor numOff sh aPtr bPtr
         RowMajor -> toBandedRowMajor numOff sh aPtr bPtr

toBandedColumnMajor ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   UnaryProxy offDiag -> size -> Ptr a -> Ptr a -> IO ()
toBandedColumnMajor numOff size aPtr bPtr = do
   let n = Shape.size size
   let k = integralFromProxy numOff
   let lda0 = k
   let lda = lda0+1
   let ldb0 = 2*k
   let ldb = ldb0+1
   copySubMatrix lda n lda aPtr ldb bPtr
   evalContT $ do
      incxPtr <- Call.cint lda0
      incyPtr <- Call.cint 1
      inczPtr <- Call.cint 0
      zPtr <- Call.number zero
      nPtr <- Call.alloca
      liftIO $ for_ (take n [0..]) $ \i -> do
         let top = i+1
         let bottom = min n (i+k+1)
         let xPtr = advancePtr aPtr ((i+1)*lda0+top+k-1)
         let yPtr = advancePtr bPtr (i*ldb0+k)
         pokeCInt nPtr (bottom-top)
         copyConjugate nPtr xPtr incxPtr (advancePtr yPtr top) incyPtr
         pokeCInt nPtr (i+k+1 - bottom)
         BlasGen.copy nPtr zPtr inczPtr (advancePtr yPtr bottom) incyPtr

toBandedRowMajor ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   UnaryProxy offDiag -> size -> Ptr a -> Ptr a -> IO ()
toBandedRowMajor numOff size aPtr bPtr = do
   let n = Shape.size size
   let k = integralFromProxy numOff
   let lda0 = k
   let lda = lda0+1
   let ldb0 = 2*k
   let ldb = ldb0+1
   copySubMatrix lda n lda aPtr ldb (advancePtr bPtr k)
   evalContT $ do
      incxPtr <- Call.cint lda0
      incyPtr <- Call.cint 1
      inczPtr <- Call.cint 0
      zPtr <- Call.number zero
      nPtr <- Call.alloca
      liftIO $ for_ (take n [0..]) $ \i -> do
         let left = max 0 (i-k)
         let xPtr = advancePtr aPtr (left*lda0+i)
         let yPtr = advancePtr bPtr (i*ldb0)
         pokeCInt nPtr (k-i+left)
         BlasGen.copy nPtr zPtr inczPtr (advancePtr yPtr i) incyPtr
         pokeCInt nPtr (i-left)
         copyConjugate nPtr xPtr incxPtr (advancePtr yPtr (left+k)) incyPtr


multiplyVector ::
   (Unary.Natural offDiag, Shape.C size, Eq size, Class.Floating a) =>
   Transposition -> BandedHermitian offDiag size a ->
   Vector size a -> Vector size a
multiplyVector transposed
   (Array (MatrixShape.BandedHermitian numOff order size) a) (Array sizeX x) =
      Array.unsafeCreateWithSize size $ \n yPtr -> do

   Call.assert "BandedHermitian.multiplyVector: shapes mismatch"
      (size == sizeX)
   let k = integralFromProxy numOff
   evalContT $ do
      let conj = transposeOrder transposed order == RowMajor
      uploPtr <- Call.char $ uploFromOrder order
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim $ k+1
      xPtr <- condConjugateToTemp conj n x
      incxPtr <- Call.cint 1
      betaPtr <- Call.number zero
      incyPtr <- Call.cint 1
      liftIO $ do
         BlasGen.hbmv uploPtr nPtr kPtr
            alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr
         when conj $ lacgv nPtr yPtr incyPtr


covariance ::
   (Shape.C size, Eq size, Class.Floating a,
    Unary.Natural sub, Unary.Natural super) =>
   Banded.Square sub super size a ->
   BandedHermitian (sub :+: super) size a
covariance a =
   case mapPair (natFromProxy,natFromProxy) $
        MatrixShape.bandedOffDiagonals $ Array.shape a of
      (sub,super) ->
         case (Proof.addNat sub super, Proof.addComm sub super) of
            (Proof.Nat, Proof.AddComm) ->
               fromUpperPart $ Banded.multiply (Banded.adjoint a) a

fromUpperPart ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   Banded.Square offDiag offDiag size a -> BandedHermitian offDiag size a
fromUpperPart (Array (MatrixShape.Banded (sub,super) order extent) a) =
   let sh = Extent.squareSize extent
       n = Shape.size sh
       kl = integralFromProxy sub
       ku = integralFromProxy super
       lda = kl+1+ku
       ldb = ku+1
   in Array.unsafeCreate (MatrixShape.BandedHermitian super order sh) $ \bPtr ->
      withForeignPtr a $ \aPtr ->
      case order of
         ColumnMajor -> copySubMatrix ldb n lda aPtr ldb bPtr
         RowMajor -> copySubMatrix ldb n lda (advancePtr aPtr kl) ldb bPtr


multiplyFull ::
   (Unary.Natural offDiag, Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Eq width, Class.Floating a) =>
   Transposition -> BandedHermitian offDiag height a ->
   Matrix.Full vert horiz height width a ->
   Matrix.Full vert horiz height width a
multiplyFull transposed a b =
   case MatrixShape.fullOrder $ Array.shape b of
      ColumnMajor -> multiplyFullSpecial transposed a b
      RowMajor -> multiplyFullGeneric transposed a b

multiplyFullSpecial ::
   (Unary.Natural offDiag, Extent.C vert, Extent.C horiz,
    Eq height, Shape.C height, Shape.C width, Class.Floating a) =>
   Transposition -> BandedHermitian offDiag height a ->
   Matrix.Full vert horiz height width a ->
   Matrix.Full vert horiz height width a
multiplyFullSpecial transposed
      (Array (MatrixShape.BandedHermitian numOff orderA sizeA) a)
      (Array (MatrixShape.Full orderB extentB) b) =
   Array.unsafeCreate (MatrixShape.Full orderB extentB) $ \cPtr -> do
      Call.assert "BandedHermitian.multiplyFull: shapes mismatch"
         (sizeA == Extent.height extentB)
      let (height,width) = Extent.dimensions extentB
      case orderB of
         ColumnMajor ->
            multiplyFullColumnMajor
               transposed numOff (height,width) orderA a b cPtr
         RowMajor ->
            multiplyFullRowMajor
               transposed numOff (height,width) orderA a b cPtr

multiplyFullColumnMajor ::
   (Unary.Natural offDiag, Shape.C height, Shape.C width, Class.Floating a) =>
   Transposition -> UnaryProxy offDiag -> (height, width) ->
   Order -> ForeignPtr a -> ForeignPtr a -> Ptr a -> IO ()
multiplyFullColumnMajor transposed numOff (height,width) order a b cPtr = do
   let n = Shape.size height
   let nrhs = Shape.size width
   let k = integralFromProxy numOff
   evalContT $ do
      uploPtr <- Call.char $ uploFromOrder order
      nPtr <- Call.cint n
      kPtr <- Call.cint k
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim $ k+1
      bPtr <- ContT $ withForeignPtr b
      incxPtr <- Call.cint 1
      betaPtr <- Call.number zero
      incyPtr <- Call.cint 1
      let pointers = take nrhs $ zip (pointerSeq n bPtr) (pointerSeq n cPtr)
      case transposeOrder transposed order of
         RowMajor -> do
            xPtr <- Call.allocaArray n
            liftIO $ for_ pointers $ \(biPtr,yPtr) -> do
               copyConjugate nPtr biPtr incxPtr xPtr incxPtr
               BlasGen.hbmv uploPtr nPtr kPtr
                  alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr
               lacgv nPtr yPtr incyPtr
         ColumnMajor ->
            liftIO $ for_ pointers $ \(xPtr,yPtr) ->
               BlasGen.hbmv uploPtr nPtr kPtr
                  alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr

multiplyFullRowMajor ::
   (Unary.Natural offDiag, Shape.C height, Shape.C width, Class.Floating a) =>
   Transposition -> UnaryProxy offDiag -> (height, width) ->
   Order -> ForeignPtr a -> ForeignPtr a -> Ptr a -> IO ()
multiplyFullRowMajor =
   error "BandedHermitian.multiplyFullRowMajor: not implemented"


multiplyFullGeneric ::
   (Unary.Natural offDiag, Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width, Eq width, Class.Floating a) =>
   Transposition -> BandedHermitian offDiag height a ->
   Matrix.Full vert horiz height width a ->
   Matrix.Full vert horiz height width a
multiplyFullGeneric transposed a b =
   let (lower,upper) = (takeStrictLower a, takeUpper a)
       (lowerT,upperT) =
         case transposed of
            Transposed -> (Banded.transpose upper, Banded.transpose lower)
            NonTransposed -> (lower,upper)
   in Banded.multiplyFull (Banded.mapExtent Extent.fromSquare lowerT) b
      `Vector.add`
      Banded.multiplyFull (Banded.mapExtent Extent.fromSquare upperT) b

takeUpper ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   BandedHermitian offDiag size a ->
   Banded.Square TypeNum.U0 offDiag size a
takeUpper =
   Array.mapShape
      (\(MatrixShape.BandedHermitian numOff order sh) ->
         MatrixShape.bandedSquare (Proxy,numOff) order sh)

takeStrictLower ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   BandedHermitian offDiag size a ->
   Banded.Square offDiag TypeNum.U0 size a
takeStrictLower (Array (MatrixShape.BandedHermitian numOff order sh) x) =
   Array.unsafeCreateWithSize
      (MatrixShape.bandedSquare (numOff,Proxy) (flipOrder order) sh) $
         \size yPtr -> evalContT $ do
   let k = integralFromProxy numOff
   nPtr <- Call.cint $ Shape.size sh
   xPtr <- ContT $ withForeignPtr x
   sizePtr <- Call.cint size
   incxPtr <- Call.cint 1
   incyPtr <- Call.cint 1
   inczPtr <- Call.cint 0
   ldbPtr <- Call.leadingDim $ k+1
   zPtr <- Call.number zero
   liftIO $ do
      copyConjugate sizePtr xPtr incxPtr yPtr incyPtr
      let offset = case order of ColumnMajor -> k; RowMajor -> 0
      BlasGen.copy nPtr zPtr inczPtr (advancePtr yPtr offset) ldbPtr


type StaticVector n = Vector (ShapeStatic.ZeroBased n)

{-
The list represents ragged rows of a sparse matrix.
-}
sumRank1 ::
   (Unary.Natural k, Shape.Indexed sh, Class.Floating a) =>
   Order -> sh ->
   [(RealOf a, (Shape.Index sh, StaticVector (Unary.Succ k) a))] ->
   BandedHermitian k sh a
sumRank1 =
   getSumRank1 $
   Class.switchFloating
      (SumRank1 $ sumRank1Aux Proxy)
      (SumRank1 $ sumRank1Aux Proxy)
      (SumRank1 $ sumRank1Aux Proxy)
      (SumRank1 $ sumRank1Aux Proxy)

newtype SumRank1 k sh a = SumRank1 {getSumRank1 :: SumRank1_ k sh (RealOf a) a}

type SumRank1_ k sh ar a =
   Order -> sh ->
   [(ar, (Shape.Index sh, StaticVector (Unary.Succ k) a))] ->
   BandedHermitian k sh a

sumRank1Aux ::
   (Unary.Natural k, Shape.Indexed sh,
    Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   UnaryProxy k -> SumRank1_ k sh ar a
sumRank1Aux numOff order size xs =
   Array.unsafeCreateWithSize
      (MatrixShape.BandedHermitian numOff order size) $
         \bSize aPtr -> evalContT $ do
   let k = integralFromProxy numOff
   let n = Shape.size size
   let lda = k+1
   uploPtr <- Call.char $ uploFromOrder order
   mPtr <- Call.cint lda
   alphaPtr <- Call.alloca
   incxPtr <- Call.cint 1
   kPtr <- Call.cint k
   ldbPtr <- Call.leadingDim k
   bSizePtr <- Call.cint bSize
   liftIO $ do
      fill zero bSize aPtr
      for_ xs $ \(alpha, (offset, Array _shX x)) ->
         withForeignPtr x $ \xPtr -> do
            let i = Shape.offset size offset
            Call.assert "BandedHermitian.sumRank1: index too large" (i+k < n)
            let bPtr = advancePtr aPtr (lda*i)
            hbr order k alpha
               uploPtr mPtr kPtr alphaPtr xPtr incxPtr bPtr incxPtr ldbPtr
      case order of
         RowMajor -> lacgv bSizePtr aPtr incxPtr
         ColumnMajor -> return ()


type HBR_ ar a =
   Order -> Int -> ar -> Ptr CChar -> Ptr CInt -> Ptr CInt ->
   Ptr a -> Ptr a -> Ptr CInt -> Ptr a -> Ptr CInt -> Ptr CInt -> IO ()

newtype HBR a = HBR {getHBR :: HBR_ (RealOf a) a}

hbr :: Class.Floating a => HBR_ (RealOf a) a
hbr = getHBR $ Class.switchFloating (HBR syr) (HBR syr) (HBR her) (HBR her)

syr :: (Class.Real a) => HBR_ a a
syr order k alpha uploPtr nPtr kPtr alphaPtr xPtr incxPtr a0Ptr incaPtr ldaPtr =
   case order of
      ColumnMajor -> do
         let aPtr = advancePtr a0Ptr k
         poke alphaPtr alpha
         BlasReal.syr uploPtr kPtr alphaPtr xPtr incxPtr aPtr ldaPtr
         poke alphaPtr . (alpha*) =<< peekElemOff xPtr k
         BlasGen.axpy nPtr alphaPtr xPtr incxPtr (advancePtr aPtr (k*k)) incaPtr
      RowMajor -> do
         let aPtr = a0Ptr
         poke alphaPtr . (alpha*) =<< peek xPtr
         BlasGen.axpy nPtr alphaPtr xPtr incxPtr aPtr incaPtr
         poke alphaPtr alpha
         BlasReal.syr uploPtr kPtr alphaPtr
            (advancePtr xPtr 1) incxPtr (advancePtr aPtr (k+1)) ldaPtr

her :: (Class.Real a) => HBR_ a (Complex a)
her order k alpha uploPtr nPtr kPtr alphaPtr xPtr incxPtr a0Ptr incaPtr ldaPtr =
   case order of
      ColumnMajor -> do
         let aPtr = advancePtr a0Ptr k
         let alphaRealPtr = castPtr alphaPtr
         poke alphaRealPtr alpha
         BlasComplex.her uploPtr kPtr alphaRealPtr xPtr incxPtr aPtr ldaPtr
         poke alphaPtr . fmap (alpha*) . conjugate =<< peekElemOff xPtr k
         BlasGen.axpy nPtr alphaPtr xPtr incxPtr (advancePtr aPtr (k*k)) incaPtr
      RowMajor -> do
         let aPtr = a0Ptr
         let alphaRealPtr = castPtr alphaPtr
         poke alphaPtr . fmap (alpha*) . conjugate =<< peek xPtr
         BlasGen.axpy nPtr alphaPtr xPtr incxPtr aPtr incaPtr
         poke alphaRealPtr alpha
         BlasComplex.her uploPtr kPtr alphaRealPtr
            (advancePtr xPtr 1) incxPtr (advancePtr aPtr (k+1)) ldaPtr
