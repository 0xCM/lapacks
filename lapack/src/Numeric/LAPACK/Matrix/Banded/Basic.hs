{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Numeric.LAPACK.Matrix.Banded.Basic (
   Banded,
   General,
   Square,
   Upper,
   Lower,
   Diagonal,
   fromList,
   squareFromList,
   lowerFromList,
   upperFromList,
   mapExtent,
   diagonal,
   takeDiagonal,
   toFull,
   toLowerTriangular,
   toUpperTriangular,
   transpose,
   adjoint,
   multiplyVector,
   multiply,
   multiplyFull,
   ) where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Matrix.Triangular.Private as TriangularPriv
import qualified Numeric.LAPACK.Matrix.Triangular.Basic as Triangular
import qualified Numeric.LAPACK.Matrix.Private as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor,ColumnMajor), transposeFromOrder, swapOnRowMajor,
          UnaryProxy, addOffDiagonals)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (zero, one)
import Numeric.LAPACK.Private
         (fill, pointerSeq, pokeCInt, copySubMatrix, copySubTrapezoid)

import qualified Numeric.BLAS.FFI.Generic as BlasGen
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
import Foreign.ForeignPtr (ForeignPtr, withForeignPtr)
import Foreign.Ptr (Ptr)
import Foreign.Storable (Storable)

import qualified Control.Monad.Trans.Maybe as MM
import qualified Control.Monad.Trans.Reader as MR
import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Monad (mzero, void)

import Data.Foldable (forM_)
import Data.Tuple.HT (swap)
import Data.Ord.HT (limit)


type Banded sub super vert horiz height width =
      Array (MatrixShape.Banded sub super vert horiz height width)

type General sub super height width =
      Array (MatrixShape.BandedGeneral sub super height width)

type Square sub super size =
      Array (MatrixShape.BandedSquare sub super size)

type Lower sub size = Square sub TypeNum.U0 size
type Upper super size = Square TypeNum.U0 super size

type Diagonal size = Square TypeNum.U0 TypeNum.U0 size


fromList ::
   (Unary.Natural sub, Unary.Natural super,
    Shape.C height, Shape.C width, Storable a) =>
   (UnaryProxy sub, UnaryProxy super) -> Order -> height -> width -> [a] ->
   General sub super height width a
fromList offDiag order height width =
   fromListGen offDiag order (Extent.general height width)

squareFromList ::
   (Unary.Natural sub, Unary.Natural super, Shape.C size, Storable a) =>
   (UnaryProxy sub, UnaryProxy super) -> Order -> size -> [a] ->
   Square sub super size a
squareFromList offDiag order size =
   fromListGen offDiag order (Extent.square size)

lowerFromList ::
   (Unary.Natural sub, Shape.C size, Storable a) =>
   UnaryProxy sub -> Order -> size -> [a] -> Lower sub size a
lowerFromList numOff order size =
   fromListGen (numOff,Proxy) order (Extent.square size)

upperFromList ::
   (Unary.Natural super, Shape.C size, Storable a) =>
   UnaryProxy super -> Order -> size -> [a] -> Upper super size a
upperFromList numOff order size =
   fromListGen (Proxy,numOff) order (Extent.square size)

fromListGen ::
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Storable a) =>
   (UnaryProxy sub, UnaryProxy super) -> Order ->
   Extent.Extent vert horiz height width -> [a] ->
   Banded sub super vert horiz height width a
fromListGen offDiag order extent =
   Array.fromList (MatrixShape.Banded offDiag order extent)


mapExtent ::
   (Extent.C vertA, Extent.C horizA) =>
   (Extent.C vertB, Extent.C horizB) =>
   Extent.Map vertA horizA vertB horizB height width ->
   Banded super sub vertA horizA height width a ->
   Banded super sub vertB horizB height width a
mapExtent f = Array.mapShape $ MatrixShape.bandedMapExtent f

transpose ::
   (Extent.C vert, Extent.C horiz) =>
   Banded sub super vert horiz height width a ->
   Banded super sub horiz vert width height a
transpose = Array.mapShape MatrixShape.bandedTranspose

adjoint ::
   (Unary.Natural super, Unary.Natural sub, Extent.C vert, Extent.C horiz,
    Shape.C width, Shape.C height, Class.Floating a) =>
   Banded sub super vert horiz height width a ->
   Banded super sub horiz vert width height a
adjoint = Vector.conjugate . transpose


diagonal :: (Shape.C sh, Class.Floating a) => Vector sh a -> Diagonal sh a
diagonal (Array sh x) =
   Array (MatrixShape.bandedSquare (Proxy,Proxy) ColumnMajor sh) x

takeDiagonal ::
   (Unary.Natural sub, Unary.Natural super, Shape.C sh, Class.Floating a) =>
   Square sub super sh a -> Vector sh a
takeDiagonal (Array (MatrixShape.Banded (sub,super) order extent) x) =
   let size = Extent.squareSize extent
       kl = integralFromProxy sub
       ku = integralFromProxy super
   in if (kl,ku) == (0,0)
        then Array size x
        else
            Array.unsafeCreateWithSize size $ \n yPtr -> evalContT $ do
               nPtr <- Call.cint n
               xPtr <- ContT $ withForeignPtr x
               let k =
                     case order of
                        RowMajor -> kl
                        ColumnMajor -> ku
               incxPtr <- Call.cint (kl+ku+1)
               incyPtr <- Call.cint 1
               liftIO $
                  BlasGen.copy nPtr (advancePtr xPtr k) incxPtr yPtr incyPtr


multiplyVector ::
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width, Eq width,
    Class.Floating a) =>
   Banded sub super vert horiz height width a ->
   Vector width a -> Vector height a
multiplyVector
   (Array (MatrixShape.Banded numOff order extent) a) (Array width x) =
      let height = Extent.height extent
      in Array.unsafeCreate height $ \yPtr -> do

   Call.assert "Banded.multiplyVector: shapes mismatch"
      (Extent.width extent == width)
   let (m,n) = MatrixShape.dimensions $ MatrixShape.Full order extent
   let (kl,ku) = MatrixShape.numOffDiagonals order numOff
   evalContT $ do
      transPtr <- Call.char $ transposeFromOrder order
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      klPtr <- Call.cint kl
      kuPtr <- Call.cint ku
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim $ kl+1+ku
      xPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint 1
      betaPtr <- Call.number zero
      incyPtr <- Call.cint 1
      liftIO $
         Private.gbmv transPtr mPtr nPtr klPtr kuPtr
            alphaPtr aPtr ldaPtr xPtr incxPtr betaPtr yPtr incyPtr


multiply ::
   (Unary.Natural subA, Unary.Natural superA,
    Unary.Natural subB, Unary.Natural superB,
    (subA :+: subB) ~ subC,
    (superA :+: superB) ~ superC,
    Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Shape.C fuse, Eq fuse,
    Class.Floating a) =>
   Banded subA superA vert horiz height fuse a ->
   Banded subB superB vert horiz fuse width a ->
   Banded subC superC vert horiz height width a
multiply
      (Array (MatrixShape.Banded numOffA orderA extentA) a)
      (Array (MatrixShape.Banded numOffB orderB extentB) b) =
   case (addOffDiagonals numOffA numOffB, Extent.fuse extentA extentB) of
      (_, Nothing) -> error "Banded.multiply: shapes mismatch"
      (((Proof.Nat, Proof.Nat), numOffC), Just extent) ->
         Array.unsafeCreate
               (MatrixShape.Banded numOffC orderB extent) $ \cPtr ->
            let (height,fuse) = Extent.dimensions extentA
                width = Extent.width extentB
            in case (orderA,orderB) of
                  (ColumnMajor,ColumnMajor) ->
                     multiplyColumnMajor ColumnMajor
                        numOffA numOffB (height,fuse,width) a b cPtr
                  (RowMajor,ColumnMajor) ->
                     multiplyColumnMajor RowMajor
                        numOffA numOffB (height,fuse,width) a b cPtr
                  (ColumnMajor,RowMajor) ->
                     multiplyColumnRowMajor
                        (swap numOffB) (swap numOffA)
                        (width,fuse,height) b a cPtr
                  (RowMajor,RowMajor) ->
                     multiplyColumnMajor ColumnMajor
                        (swap numOffB) (swap numOffA)
                        (width,fuse,height) b a cPtr

multiplyColumnMajor ::
   (Unary.Natural subA, Unary.Natural superA,
    Unary.Natural subB, Unary.Natural superB,
    Shape.C height, Shape.C width, Shape.C fuse,
    Class.Floating a) =>
   Order ->
   (UnaryProxy subA, UnaryProxy superA) ->
   (UnaryProxy subB, UnaryProxy superB) ->
   (height, fuse, width) ->
   ForeignPtr a -> ForeignPtr a -> Ptr a -> IO ()
multiplyColumnMajor orderA (subA,superA) (subB,superB)
      (height,fuse,width) a b cPtr = do
   let m = Shape.size height
   let k = Shape.size fuse
   let n = Shape.size width
   let (kla,kua) = (integralFromProxy subA, integralFromProxy superA)
   let (klb,kub) = (integralFromProxy subB, integralFromProxy superB)
   let ku = kua+kub
   let kl = kla+klb
   let lda0 = kla+kua
   let ldb0 = klb+kub
   let ldc0 = lda0+ldb0
   let lda = lda0+1
   let ldc = ldc0+1
   evalContT $ do
      transPtr <- Call.char $ transposeFromOrder orderA
      mPtr <- Call.alloca
      nPtr <- Call.alloca
      klPtr <- Call.alloca
      kuPtr <- Call.alloca
      let ((miPtr,kliPtr),(niPtr,kuiPtr)) =
            swapOnRowMajor orderA ((mPtr,klPtr),(nPtr,kuPtr))
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim lda
      bPtr <- ContT $ withForeignPtr b
      incxPtr <- Call.cint 1
      betaPtr <- Call.number zero
      incyPtr <- Call.cint 1
      liftIO $
         forM_ (take n [0..]) $ \i -> do
            let top = max 0 (i-ku)
            let bottom = min m (i+kl+1)
            let left = max 0 (i-kub)
            let right = min k (i+klb+1)
            pokeCInt miPtr $ max 0 $ bottom-top
            pokeCInt niPtr $ max 0 $ right-left
            let d = top-left; kli = kla-d; kui = kua+d
            pokeCInt kuiPtr kui
            pokeCInt kliPtr kli
            let j0 = i*ldc
            let j1 = i*ldc0 + top+ku
            let j2 = i*ldc0 + bottom+ku
            fill zero (j1-j0) (advancePtr cPtr j0)
            let aOffset =
                  case orderA of
                     ColumnMajor -> left
                     RowMajor -> top
            Private.gbmv transPtr mPtr nPtr klPtr kuPtr
               alphaPtr
               (advancePtr aPtr (aOffset*lda)) ldaPtr
               (advancePtr bPtr (i*ldb0 + left+kub)) incxPtr
               betaPtr
               (advancePtr cPtr j1) incyPtr
            fill zero (j0+ldc-j2) (advancePtr cPtr j2)

multiplyColumnRowMajor ::
   (Unary.Natural subA, Unary.Natural superA,
    Unary.Natural subB, Unary.Natural superB,
    Shape.C height, Shape.C width, Shape.C fuse,
    Class.Floating a) =>
   (UnaryProxy subA, UnaryProxy superA) ->
   (UnaryProxy subB, UnaryProxy superB) ->
   (height, fuse, width) ->
   ForeignPtr a -> ForeignPtr a -> Ptr a -> IO ()
multiplyColumnRowMajor (subA,superA) (subB,superB)
      (height,fuse,width) a b cPtr = do
   let m = Shape.size height
   let k = Shape.size fuse
   let n = Shape.size width
   let (kla,kua) = (integralFromProxy subA, integralFromProxy superA)
   let (klb,kub) = (integralFromProxy subB, integralFromProxy superB)
   let ku = kua+kub
   let kl = kla+klb
   let lda0 = kla+kua
   let ldb0 = klb+kub
   let ldc0 = kl+ku
   let ldc = ldc0+1
   fill zero (ldc*n) cPtr
   evalContT $ do
      mPtr <- Call.alloca
      nPtr <- Call.alloca
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      bPtr <- ContT $ withForeignPtr b
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 1
      ldc0Ptr <- Call.leadingDim $ ldc0 + if ldb0==0 then 1 else 0
      liftIO $
         forM_ (take k [0..]) $ \i -> do
            let top = max 0 (i-kua)
            let bottom = min m (i+kla+1)
            let left = max 0 (i-klb)
            let right = min n (i+kub+1)
            pokeCInt mPtr $ max 0 $ bottom-top
            pokeCInt nPtr $ max 0 $ right-left
            BlasGen.geru mPtr nPtr alphaPtr
               (advancePtr aPtr (i*lda0+top+kua)) incxPtr
               (advancePtr bPtr (i*ldb0+left+klb)) incyPtr
               (advancePtr cPtr (left*ldc0+top+ku)) ldc0Ptr


multiplyFull ::
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Shape.C fuse, Eq fuse,
    Class.Floating a) =>
   Banded sub super vert horiz height fuse a ->
   Matrix.Full vert horiz fuse width a -> Matrix.Full vert horiz height width a
multiplyFull
      (Array (MatrixShape.Banded numOff orderA extentA) a)
      (Array (MatrixShape.Full orderB extentB) b) =
   case Extent.fuse extentA extentB of
      Nothing -> error "Banded.multiplyFull: shapes mismatch"
      Just extent ->
         Array.unsafeCreate (MatrixShape.Full orderB extent) $ \cPtr ->
            let (height,fuse) = Extent.dimensions extentA
                width = Extent.width extentB
            in case orderB of
                  ColumnMajor ->
                     multiplyFullColumnMajor
                        numOff (height,fuse,width) orderA extentA a b cPtr
                  RowMajor ->
                     multiplyFullRowMajor
                        numOff (height,fuse,width) orderA a b cPtr

multiplyFullColumnMajor ::
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz,
    Shape.C height, Shape.C width, Shape.C fuse,
    Class.Floating a) =>
   (UnaryProxy sub, UnaryProxy super) ->
   (height, fuse, width) ->
   Order -> Extent.Extent vert horiz height fuse ->
   ForeignPtr a -> ForeignPtr a -> Ptr a -> IO ()
multiplyFullColumnMajor numOff (height,fuse,width) orderA extentA a b cPtr = do
   let (m,n) = MatrixShape.dimensions $ MatrixShape.Full orderA extentA
   let k = Shape.size width
   let (kl,ku) = MatrixShape.numOffDiagonals orderA numOff
   evalContT $ do
      transPtr <- Call.char $ transposeFromOrder orderA
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      klPtr <- Call.cint kl
      kuPtr <- Call.cint ku
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim $ kl+1+ku
      bPtr <- ContT $ withForeignPtr b
      incxPtr <- Call.cint 1
      betaPtr <- Call.number zero
      incyPtr <- Call.cint 1
      liftIO $
         forM_ (take k $
                zip (pointerSeq (Shape.size fuse) bPtr)
                    (pointerSeq (Shape.size height) cPtr)) $
            \(xPtr,yPtr) ->
               Private.gbmv transPtr mPtr nPtr klPtr kuPtr
                  alphaPtr aPtr ldaPtr xPtr incxPtr
                  betaPtr yPtr incyPtr

multiplyFullRowMajor ::
   (Unary.Natural sub, Unary.Natural super,
    Shape.C height, Shape.C width, Shape.C fuse,
    Class.Floating a) =>
   (UnaryProxy sub, UnaryProxy super) ->
   (height, fuse, width) ->
   Order -> ForeignPtr a -> ForeignPtr a -> Ptr a -> IO ()
multiplyFullRowMajor (sub,super) (height,fuse,width) orderA a b cPtr = do
   let m = Shape.size height
   let n = Shape.size fuse
   let k = Shape.size width
   let kl = integralFromProxy sub
   let ku = integralFromProxy super
   let lda0 = kl+ku
   let lda = lda0+1
   evalContT $ do
      transPtr <- Call.char 'N'
      kPtr <- Call.cint k
      dPtr <- Call.alloca
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      bPtr <- ContT $ withForeignPtr b
      ldbPtr <- Call.leadingDim k
      incxPtr <- Call.cint $
         case orderA of
            RowMajor -> 1
            ColumnMajor -> max 1 lda0
      betaPtr <- Call.number zero
      incyPtr <- Call.cint 1
      liftIO $
         forM_ (take m $ zip [0..] $
                zip (pointerSeq lda aPtr) (pointerSeq k cPtr)) $
            \(i,(xPtr,yPtr)) -> do
               let firstRow = limit (0,n) (i-kl)
               let last1Row = limit (0,n) (i+ku+1)
               let biPtr = advancePtr bPtr (firstRow*k)
               let xOffset =
                     case orderA of
                        RowMajor -> firstRow-i+kl
                        ColumnMajor -> (firstRow-i)*lda0+ku
               let xiPtr = advancePtr xPtr xOffset
               pokeCInt dPtr $ last1Row - firstRow
               Private.gemv transPtr kPtr dPtr
                  alphaPtr biPtr ldbPtr xiPtr incxPtr
                  betaPtr yPtr incyPtr


toLowerTriangular ::
   (Unary.Natural sub, Shape.C sh, Class.Floating a) =>
   Lower sub sh a -> Triangular.Lower sh a
toLowerTriangular =
   Triangular.transpose . toUpperTriangular . transpose

toUpperTriangular ::
   (Unary.Natural super, Shape.C sh, Class.Floating a) =>
   Upper super sh a -> Triangular.Upper sh a
toUpperTriangular (Array (MatrixShape.Banded (_sub,super) order extent) a) =
   let size = Extent.squareSize extent
   in Array.unsafeCreateWithSize
         (MatrixShape.Triangular MatrixShape.NonUnit MatrixShape.upper
            order size) $
      TriangularPriv.fromBanded
         (integralFromProxy super) order (Shape.size size) a

toFull ::
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Banded sub super vert horiz height width a ->
   Matrix.Full vert horiz height width a
toFull (Array (MatrixShape.Banded (sub,super) order extent) a) =
   Array.unsafeCreateWithSize (MatrixShape.Full order extent) $ \bSize bPtr ->
   withForeignPtr a $ \aPtr -> do
      let (height,width) = Extent.dimensions extent
      fill zero bSize bPtr
      case order of
         ColumnMajor -> toFullColumnMajor (sub,super) (height,width) aPtr bPtr
         RowMajor -> toFullColumnMajor (super,sub) (width,height) aPtr bPtr

toFullColumnMajor ::
   (Unary.Natural sub, Unary.Natural super, Shape.C height, Shape.C width,
    Class.Floating a) =>
   (UnaryProxy sub, UnaryProxy super) -> (height,width) ->
   Ptr a -> Ptr a -> IO ()
toFullColumnMajor (sub,super) (height,width) aPtr bPtr = do
   let m = Shape.size height
   let n = Shape.size width
   let kl = integralFromProxy sub
   let ku = integralFromProxy super
   let lda0 = kl+ku
   let lda = lda0+1

   void $ MM.runMaybeT $ flip MR.runReaderT n $
      if m > lda0
         then do -- diagonal stripe
            let col0 = ku
            withRightBound col0 $ \col ->
               copyUpperTrapezoid (col+kl) col lda0 (advancePtr aPtr ku) m bPtr
            let col1 = m-kl
            withRightBound col1 $ \col ->
               copySubMatrix lda (col-col0)
                  lda (advancePtr aPtr (col0*lda))
                  (m+1) (advancePtr bPtr (col0*m))
            let col2 = m+ku
            withRightBound col2 $ \col ->
               copySubTrapezoid 'L' lda0 (col-col1)
                  lda0 (advancePtr aPtr (col1*lda))
                  m (advancePtr bPtr (col1*m+m-lda0))
         else do -- full block in the middle
            let col0 = max 0 $ m-kl
            withRightBound col0 $ \col ->
               copyUpperTrapezoid (col+kl) col lda0 (advancePtr aPtr ku) m bPtr
            let col1 = ku
            withRightBound col1 $ \col ->
               copySubMatrix m (col-col0)
                  lda0 (advancePtr aPtr (col0*lda+(col1-col0)))
                  m (advancePtr bPtr (col0*m))
            let col2 = m+ku
            withRightBound col2 $ \col ->
               copySubTrapezoid 'L' m (col-col1)
                  lda0 (advancePtr aPtr (ku*lda))
                  m (advancePtr bPtr (ku*m))

withRightBound ::
   Int -> (Int -> IO a) -> MR.ReaderT Int (MM.MaybeT IO) a
withRightBound col act = do
   n <- MR.ask
   if n<=col
     then liftIO (act n) >> mzero
     else liftIO (act col)

copyUpperTrapezoid ::
   (Class.Floating a) =>
   Int -> Int -> Int -> Ptr a -> Int -> Ptr a -> IO ()
copyUpperTrapezoid m n lda aPtr ldb bPtr = do
   let d = m-n
   copySubMatrix d n lda aPtr ldb bPtr
   copySubTrapezoid 'U' n n
      lda (advancePtr aPtr d)
      ldb (advancePtr bPtr d)
