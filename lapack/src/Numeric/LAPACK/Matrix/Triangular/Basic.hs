{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
module Numeric.LAPACK.Matrix.Triangular.Basic (
   Triangular, MatrixShape.UpLo,
   Upper, UnitUpper,
   Lower, UnitLower,
   Symmetric, Diagonal,
   fromList, autoFromList,
   lowerFromList, autoLowerFromList,
   upperFromList, autoUpperFromList,
   symmetricFromList, autoSymmetricFromList,
   diagonalFromList, autoDiagonalFromList,
   relaxUnitDiagonal, strictNonUnitDiagonal,
   asDiagonal, asLower, asUpper, asSymmetric,
   forceUnitDiagonal, forceNonUnitDiagonal,
   identity,
   diagonal,
   takeDiagonal,
   transpose,
   adjoint,

   toSquare,
   takeUpper,
   takeLower,

   fromLowerRowMajor, toLowerRowMajor,
   fromUpperRowMajor, toUpperRowMajor,
   forceOrder, adaptOrder,

   add, sub,

   Tri.PowerDiag,
   multiplyVector,
   square, squareGeneric,
   multiply,
   multiplyFull,
   ) where

import qualified Numeric.LAPACK.Matrix.Symmetric.Private as Symmetric
import qualified Numeric.LAPACK.Matrix.Triangular.Private as Tri
import qualified Numeric.LAPACK.Matrix.Basic as Basic
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix.Triangular.Private
         (Triangular, FlexDiagonal, diagonalPointers, diagonalPointerPairs,
          pack, packRect, unpack, unpackZero, unpackToTemp)
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor,ColumnMajor),
          flipOrder, transposeFromOrder, uploFromOrder, uploOrder,
          Unit(Unit), NonUnit(NonUnit), charFromTriDiag)
import Numeric.LAPACK.Matrix.Private
         (Full, Square, ZeroInt, zeroInt, Conjugation(NonConjugated))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (zero, one)
import Numeric.LAPACK.Private (fill, copyBlock)

import qualified Numeric.LAPACK.FFI.Complex as LapackComplex
import qualified Numeric.BLAS.FFI.Real as BlasReal
import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.C.Types (CChar, CInt)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Ptr (Ptr)
import Foreign.Storable (Storable, poke, peek)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)

import Data.Foldable (forM_)


type Lower sh = FlexLower NonUnit sh
type Upper sh = FlexUpper NonUnit sh
type Symmetric sh = Array (MatrixShape.Symmetric sh)
type Diagonal sh = FlexDiagonal NonUnit sh

type FlexLower diag sh = Array (MatrixShape.LowerTriangular diag sh)
type FlexUpper diag sh = Array (MatrixShape.UpperTriangular diag sh)
type FlexSymmetric diag sh = Array (MatrixShape.FlexSymmetric diag sh)

type UnitLower sh = Array (MatrixShape.LowerTriangular Unit sh)
type UnitUpper sh = Array (MatrixShape.UpperTriangular Unit sh)

transpose ::
   (MatrixShape.Content lo, MatrixShape.Content up,
    MatrixShape.TriDiag diag) =>
   Triangular lo diag up sh a -> Triangular up diag lo sh a
transpose (Array sh a) =
   Array (MatrixShape.triangularTranspose sh) a

adjoint ::
   (MatrixShape.Content lo, MatrixShape.Content up,
    MatrixShape.TriDiag diag, Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Triangular up diag lo sh a
adjoint = Vector.conjugate . transpose


fromList ::
   (MatrixShape.Content lo, MatrixShape.Content up, Shape.C sh, Storable a) =>
   Order -> sh -> [a] -> Triangular lo NonUnit up sh a
fromList order sh =
   Array.fromList (MatrixShape.Triangular NonUnit MatrixShape.autoUplo order sh)

lowerFromList :: (Shape.C sh, Storable a) => Order -> sh -> [a] -> Lower sh a
lowerFromList = fromList

upperFromList :: (Shape.C sh, Storable a) => Order -> sh -> [a] -> Upper sh a
upperFromList = fromList

symmetricFromList ::
   (Shape.C sh, Storable a) => Order -> sh -> [a] -> Symmetric sh a
symmetricFromList = fromList

diagonalFromList ::
   (Shape.C sh, Storable a) => Order -> sh -> [a] -> Diagonal sh a
diagonalFromList = fromList


autoFromList ::
   (MatrixShape.Content lo, MatrixShape.Content up, Storable a) =>
   Order -> [a] -> Triangular lo NonUnit up ZeroInt a
autoFromList order xs =
   let n = length xs
       triSize = MatrixShape.triangleExtent "Triangular.autoFromList" n
       uplo = MatrixShape.autoUplo
       size = MatrixShape.caseDiagUpLoSym uplo n triSize triSize triSize
   in Array.fromList
         (MatrixShape.Triangular
            MatrixShape.autoDiag uplo order (zeroInt size))
         xs

autoLowerFromList :: (Storable a) => Order -> [a] -> Lower ZeroInt a
autoLowerFromList = autoFromList

autoUpperFromList :: (Storable a) => Order -> [a] -> Upper ZeroInt a
autoUpperFromList = autoFromList

autoSymmetricFromList :: (Storable a) => Order -> [a] -> Symmetric ZeroInt a
autoSymmetricFromList = autoFromList

autoDiagonalFromList :: (Storable a) => Order -> [a] -> Diagonal ZeroInt a
autoDiagonalFromList = autoFromList


asDiagonal :: FlexDiagonal diag sh a -> FlexDiagonal diag sh a
asDiagonal = id

asLower :: FlexLower diag sh a -> FlexLower diag sh a
asLower = id

asUpper :: FlexUpper diag sh a -> FlexUpper diag sh a
asUpper = id

asSymmetric :: FlexSymmetric diag sh a -> FlexSymmetric diag sh a
asSymmetric = id

forceUnitDiagonal :: Triangular lo Unit up sh a -> Triangular lo Unit up sh a
forceUnitDiagonal = id

forceNonUnitDiagonal ::
   Triangular lo NonUnit up sh a -> Triangular lo NonUnit up sh a
forceNonUnitDiagonal = id


toSquare ::
   (MatrixShape.Content lo, MatrixShape.Content up,
    Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Square sh a
toSquare (Array (MatrixShape.Triangular _diag uplo order sh) a) =
   Array.unsafeCreateWithSize (MatrixShape.square order sh) $ \size bPtr ->
      let n = Shape.size sh
      in withForeignPtr a $ \aPtr ->
            MatrixShape.caseDiagUpLoSym uplo
               (do
                  fill zero size bPtr
                  evalContT $ do
                     nPtr <- Call.cint n
                     incxPtr <- Call.cint 1
                     incyPtr <- Call.cint (n+1)
                     liftIO $ BlasGen.copy nPtr aPtr incxPtr bPtr incyPtr)
               (unpackZero order n aPtr bPtr)
               (unpackZero (flipOrder order) n aPtr bPtr)
               (Symmetric.unpack NonConjugated order n aPtr bPtr)

takeLower ::
   (Extent.C horiz, Shape.C height, Shape.C width, Class.Floating a) =>
   Full Extent.Small horiz height width a -> Lower height a
takeLower =
   Tri.takeLower (MatrixShape.NonUnit, const $ const $ const $ return ())

takeUpper ::
   (Extent.C vert, Shape.C height, Shape.C width, Class.Floating a) =>
   Full vert Extent.Small height width a -> Upper width a
takeUpper (Array (MatrixShape.Full order extent) a) =
   let (height,width) = Extent.dimensions extent
       m = Shape.size height
       n = Shape.size width
       k = case order of RowMajor -> n; ColumnMajor -> m
   in Array.unsafeCreate
         (MatrixShape.Triangular MatrixShape.NonUnit
            MatrixShape.upper order width) $ \bPtr ->
      withForeignPtr a $ \aPtr -> packRect order n k aPtr bPtr


fromLowerRowMajor ::
   (Shape.C sh, Class.Floating a) =>
   Array (Shape.Triangular Shape.Lower sh) a -> Lower sh a
fromLowerRowMajor =
   Array.mapShape
      (MatrixShape.Triangular MatrixShape.NonUnit MatrixShape.lower RowMajor .
       Shape.triangularSize)

fromUpperRowMajor ::
   (Shape.C sh, Class.Floating a) =>
   Array (Shape.Triangular Shape.Upper sh) a -> Upper sh a
fromUpperRowMajor =
   Array.mapShape
      (MatrixShape.Triangular MatrixShape.NonUnit MatrixShape.upper RowMajor .
       Shape.triangularSize)

toLowerRowMajor ::
   (Shape.C sh, Class.Floating a) =>
   Lower sh a -> Array (Shape.Triangular Shape.Lower sh) a
toLowerRowMajor =
   Array.mapShape (Shape.Triangular Shape.Lower . MatrixShape.triangularSize)
   .
   forceOrder MatrixShape.RowMajor

toUpperRowMajor ::
   (Shape.C sh, Class.Floating a) =>
   Upper sh a -> Array (Shape.Triangular Shape.Upper sh) a
toUpperRowMajor =
   Array.mapShape (Shape.Triangular Shape.Upper . MatrixShape.triangularSize)
   .
   forceOrder MatrixShape.RowMajor

{-
This is not maximally efficient.
It fills up a whole square.
This wastes memory but enables more regular memory access patterns.
Additionally, it fills unused parts of the square with zero or mirrored values.
-}
forceOrder ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Shape.C sh, Class.Floating a) =>
   Order -> Triangular lo diag up sh a -> Triangular lo diag up sh a
forceOrder newOrder =
   Tri.getMap $
   MatrixShape.switchDiagUpLoSym
      (Tri.Map $
       Array.mapShape (\sh -> sh{MatrixShape.triangularOrder = newOrder}))
      (forceOrderMap newOrder takeUpper)
      (forceOrderMap newOrder takeLower)
      (forceOrderMap newOrder $
         Array.mapShape
            (\sh -> sh{MatrixShape.triangularUplo = MatrixShape.autoUplo})
         .
         takeUpper)

forceOrderMap ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Shape.C sh, Class.Floating a) =>
   Order ->
   (Square sh a -> Triangular lo NonUnit up sh a) ->
   Tri.Map diag sh a lo up
forceOrderMap newOrder f = Tri.Map $ \a ->
   if MatrixShape.triangularOrder (Array.shape a) == newOrder
      then a
      else uncheckedRelaxNonUnitDiagonal $
           f $ Basic.forceOrder newOrder $ toSquare a

uncheckedRelaxNonUnitDiagonal ::
   (MatrixShape.TriDiag diag) =>
   Triangular lo NonUnit up sh a -> Triangular lo diag up sh a
uncheckedRelaxNonUnitDiagonal =
   Array.mapShape (\sh -> sh{MatrixShape.triangularDiag = MatrixShape.autoDiag})

{- |
@adaptOrder x y@ contains the data of @y@ with the layout of @x@.
-}
adaptOrder ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a ->
   Triangular lo diag up sh a ->
   Triangular lo diag up sh a
adaptOrder x = forceOrder (MatrixShape.triangularOrder $ Array.shape x)

add, sub ::
   (MatrixShape.Content lo, MatrixShape.Content up,
    Eq lo, Eq up, Eq sh, Shape.C sh, Class.Floating a) =>
   Triangular lo NonUnit up sh a ->
   Triangular lo NonUnit up sh a ->
   Triangular lo NonUnit up sh a
add x y = Vector.add (adaptOrder y x) y
sub x y = Vector.sub (adaptOrder y x) y


identity ::
   (MatrixShape.Content lo, MatrixShape.Content up,
    Shape.C sh, Class.Floating a) =>
   Order -> sh -> Triangular lo Unit up sh a
identity order sh =
   let (realOrder, uplo) = autoUploOrder order
   in Array.unsafeCreateWithSize (MatrixShape.Triangular Unit uplo order sh) $
         \size aPtr -> do
      let n = Shape.size sh
      let fillTriangle = do
            fill zero size aPtr
            mapM_ (flip poke one) (diagonalPointers realOrder n aPtr)
      MatrixShape.caseDiagUpLoSym uplo
         (fill one n aPtr)
         fillTriangle
         fillTriangle
         fillTriangle

diagonal, diagonalAux ::
   (MatrixShape.Content lo, MatrixShape.Content up,
    Shape.C sh, Class.Floating a) =>
   Order -> Vector sh a -> Triangular lo NonUnit up sh a
diagonal order x@(Array sh xPtr) =
   let uplo = MatrixShape.autoUplo
   in MatrixShape.caseDiagUpLoSym uplo
         (Array (MatrixShape.Triangular NonUnit uplo order sh) xPtr)
         (diagonalAux order x)
         (diagonalAux order x)
         (diagonalAux order x)

diagonalAux order (Array sh x) =
   let (realOrder, uplo) = autoUploOrder order
   in Array.unsafeCreateWithSize
         (MatrixShape.Triangular NonUnit uplo order sh) $
            \size aPtr -> do
      let n = Shape.size sh
      fill zero size aPtr
      withForeignPtr x $ \xPtr ->
         forM_ (diagonalPointerPairs realOrder n xPtr aPtr) $
            \(srcPtr,dstPtr) -> poke dstPtr =<< peek srcPtr


takeDiagonal, takeDiagonalAux ::
   (MatrixShape.Content lo, MatrixShape.Content up,
    Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Vector sh a
takeDiagonal a@(Array (MatrixShape.Triangular _diag uplo _order sh) aPtr) =
   MatrixShape.caseDiagUpLoSym uplo
      (Array sh aPtr)
      (takeDiagonalAux a)
      (takeDiagonalAux a)
      (takeDiagonalAux a)

takeDiagonalAux (Array (MatrixShape.Triangular _diag uplo order sh) a) =
   Array.unsafeCreate sh $ \xPtr ->
   withForeignPtr a $ \aPtr ->
      mapM_
         (\(dstPtr,srcPtr) -> poke dstPtr =<< peek srcPtr)
         (diagonalPointerPairs (uploOrder uplo order) (Shape.size sh) xPtr aPtr)

relaxUnitDiagonal ::
   (MatrixShape.TriDiag diag) =>
   Triangular lo Unit up sh a -> Triangular lo diag up sh a
relaxUnitDiagonal = Array.mapShape MatrixShape.relaxUnitDiagonal

strictNonUnitDiagonal ::
   (MatrixShape.TriDiag diag) =>
   Triangular lo diag up sh a -> Triangular lo NonUnit up sh a
strictNonUnitDiagonal = Array.mapShape MatrixShape.strictNonUnitDiagonal


multiplyVector ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Shape.C sh, Eq sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Vector sh a -> Vector sh a
multiplyVector =
   Tri.getMultiplyRight $
   MatrixShape.switchDiagUpLoSym
      (Tri.MultiplyRight $
       Tri.multiplyDiagonal
         "multiplyVector.diagonal: sizes mismatch"
         Array.shape
         (Vector.mul . takeDiagonal))
      (Tri.MultiplyRight multiplyVectorTriangular)
      (Tri.MultiplyRight multiplyVectorTriangular)
      (Tri.MultiplyRight multiplyVectorTriangular)

multiplyVectorTriangular ::
   (MatrixShape.UpLoSym lo up, MatrixShape.TriDiag diag,
    Shape.C sh, Eq sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Vector sh a -> Vector sh a
multiplyVectorTriangular
   (Array (MatrixShape.Triangular diag uplo order shA) a) (Array shX x) =
      Array.unsafeCreate shX $ \yPtr -> do
   Call.assert "Triangular.multiplyVector: width shapes mismatch" (shA == shX)
   let n = Shape.size shA
   evalContT $ do
      uploPtr <- Call.char $ uploFromOrder $ uploOrder uplo order
      transPtr <- Call.char $ transposeFromOrder order
      diagPtr <- Call.char $ charFromTriDiag diag
      nPtr <- Call.cint n
      aPtr <- ContT $ withForeignPtr a
      xPtr <- ContT $ withForeignPtr x
      alphaPtr <- Call.number one
      betaPtr <- Call.number zero
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 1
      let runTPMV = do
            copyBlock n xPtr yPtr
            BlasGen.tpmv uploPtr transPtr diagPtr nPtr aPtr yPtr incyPtr
      liftIO $
         MatrixShape.caseUpLoSym uplo
            runTPMV
            runTPMV
            (spmv uploPtr nPtr alphaPtr aPtr xPtr incxPtr betaPtr yPtr incyPtr)


newtype SPMV a =
   SPMV {
      getSPMV ::
         Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a ->
         Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
   }

spmv :: Class.Floating a =>
   Ptr CChar -> Ptr CInt -> Ptr a -> Ptr a ->
   Ptr a -> Ptr CInt -> Ptr a -> Ptr a -> Ptr CInt -> IO ()
spmv =
   getSPMV $
   Class.switchFloating
      (SPMV BlasReal.spmv) (SPMV BlasReal.spmv)
      (SPMV LapackComplex.spmv) (SPMV LapackComplex.spmv)


square ::
   (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
    Shape.C sh, Eq sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Triangular lo diag up sh a
square =
   Tri.getMap $
   MatrixShape.switchDiagUpLo
      (Tri.Map squareDiagonal)
      (Tri.Map squareTriangular)
      (Tri.Map squareTriangular)


{- |
Include symmetric matrices.
However, symmetric matrices do not preserve unit diagonals.
-}
squareGeneric ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Shape.C sh, Eq sh, Class.Floating a) =>
   Triangular lo diag up sh a ->
   Triangular lo (Tri.PowerDiag lo up diag) up sh a
squareGeneric =
   Tri.getPower $
   MatrixShape.switchDiagUpLoSym
      (Tri.Power squareDiagonal)
      (Tri.Power squareTriangular)
      (Tri.Power squareTriangular)
      (Tri.Power $ squareSymmetric . strictNonUnitDiagonal)


squareDiagonal ::
   (MatrixShape.TriDiag diag,
    Shape.C sh, Eq sh, Class.Floating a) =>
   FlexDiagonal diag sh a -> FlexDiagonal diag sh a
squareDiagonal =
   getMapDiag $
   MatrixShape.switchTriDiag (MapDiag id) (MapDiag $ \a -> Vector.mul a a)

newtype MapDiag lo up sh a diag =
   MapDiag {
      getMapDiag ::
         Triangular lo diag up sh a ->
         Triangular lo diag up sh a
   }

squareTriangular ::
   (MatrixShape.UpLo lo up, MatrixShape.TriDiag diag,
    Shape.C sh, Eq sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Triangular lo diag up sh a
squareTriangular
   (Array shape@(MatrixShape.Triangular diag uplo order sh) a) =
      Array.unsafeCreate shape $ \bpPtr -> do
   let n = Shape.size sh
   evalContT $ do
      sidePtr <- Call.char 'L'
      let realOrder = uploOrder uplo order
      uploPtr <- Call.char $ uploFromOrder realOrder
      transPtr <- Call.char 'N'
      diagPtr <- Call.char $ charFromTriDiag diag
      nPtr <- Call.cint n
      ldPtr <- Call.leadingDim n
      aPtr <- unpackToTemp (unpack realOrder) n a
      bPtr <- unpackToTemp (unpackZero realOrder) n a
      alphaPtr <- Call.number one
      liftIO $ do
         BlasGen.trmm sidePtr uploPtr transPtr diagPtr
            nPtr nPtr alphaPtr aPtr ldPtr bPtr ldPtr
         pack realOrder n bPtr bpPtr

squareSymmetric ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Symmetric sh a -> Symmetric sh a
squareSymmetric (Array shape@(MatrixShape.Triangular _diag _uplo order sh) a) =
   Array.unsafeCreate shape $
      Symmetric.square NonConjugated order (Shape.size sh) a


multiply ::
   (MatrixShape.DiagUpLo lo up, MatrixShape.TriDiag diag,
    Shape.C sh, Eq sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Triangular lo diag up sh a ->
   Triangular lo diag up sh a
multiply =
   getMultiply $
   MatrixShape.switchDiagUpLo
      (Multiply $
       Tri.multiplyDiagonal
         "multiply.diagonal: sizes mismatch"
         (MatrixShape.triangularSize . Array.shape)
         (\a b ->
            Array.mapShape
               (MatrixShape.Triangular
                  MatrixShape.autoDiag MatrixShape.autoUplo
                  (MatrixShape.triangularOrder $ Array.shape b)) $
            Vector.mul (takeDiagonal a) (takeDiagonal b)))
      (Multiply multiplyTriangular)
      (Multiply multiplyTriangular)

newtype Multiply diag sh a lo up =
   Multiply {
      getMultiply ::
         Triangular lo diag up sh a ->
         Triangular lo diag up sh a -> Triangular lo diag up sh a
   }

multiplyTriangular ::
   (MatrixShape.UpLo lo up, MatrixShape.TriDiag diag,
    Shape.C sh, Eq sh, Class.Floating a) =>
   Triangular lo diag up sh a ->
   Triangular lo diag up sh a -> Triangular lo diag up sh a
multiplyTriangular
   (Array        (MatrixShape.Triangular diag uploA orderA shA) a)
   (Array shapeB@(MatrixShape.Triangular _diag uploB orderB shB) b) =
      Array.unsafeCreate shapeB $ \cpPtr -> do
   Call.assert "Triangular.multiply: width shapes mismatch" (shA == shB)
   let n = Shape.size shA
   evalContT $ do
      let (side,trans) =
            case orderB of
               ColumnMajor -> ('L', orderA)
               RowMajor -> ('R', flipOrder orderA)
      sidePtr <- Call.char side
      let realOrderA = uploOrder uploA orderA
      let realOrderB = uploOrder uploB orderB
      uploPtr <- Call.char $ uploFromOrder realOrderA
      transPtr <- Call.char $ transposeFromOrder trans
      diagPtr <- Call.char $ charFromTriDiag diag
      nPtr <- Call.cint n
      ldPtr <- Call.leadingDim n
      aPtr <- unpackToTemp (unpack realOrderA) n a
      bPtr <- unpackToTemp (unpackZero realOrderB) n b
      alphaPtr <- Call.number one
      liftIO $ do
         BlasGen.trmm sidePtr uploPtr transPtr diagPtr
            nPtr nPtr alphaPtr aPtr ldPtr bPtr ldPtr
         pack realOrderB n bPtr cpPtr


multiplyFull ::
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width,
    Class.Floating a) =>
   Triangular lo diag up height a ->
   Full vert horiz height width a ->
   Full vert horiz height width a
multiplyFull =
   Tri.getMultiplyRight $
   MatrixShape.switchDiagUpLoSym
      (Tri.MultiplyRight $
       Tri.multiplyDiagonal
         "multiplyFull.diagonal: sizes mismatch"
         (MatrixShape.fullHeight . Array.shape)
         (Basic.scaleRows . takeDiagonal))
      (Tri.MultiplyRight multiplyFullTriangular)
      (Tri.MultiplyRight multiplyFullTriangular)
      (Tri.MultiplyRight multiplyFullTriangular)

multiplyFullTriangular ::
   (MatrixShape.UpLoSym lo up, MatrixShape.TriDiag diag,
    Extent.C vert, Extent.C horiz,
    Shape.C height, Eq height, Shape.C width,
    Class.Floating a) =>
   Triangular lo diag up height a ->
   Full vert horiz height width a ->
   Full vert horiz height width a
multiplyFullTriangular
   (Array        (MatrixShape.Triangular diag uploA orderA shA) a)
   (Array shapeB@(MatrixShape.Full orderB extentB) b) =
      Array.unsafeCreateWithSize shapeB $ \size cPtr -> do
   let (height,width) = Extent.dimensions extentB
   Call.assert "Triangular.multiplyFull: shapes mismatch" (shA == height)
   let m0 = Shape.size height
   let n0 = Shape.size width
   evalContT $ do
      let (side,trans,(m,n)) =
            case orderB of
               ColumnMajor -> ('L', orderA, (m0,n0))
               RowMajor -> ('R', flipOrder orderA, (n0,m0))
      sidePtr <- Call.char side
      let realOrderA = uploOrder uploA orderA
      uploPtr <- Call.char $ uploFromOrder realOrderA
      transPtr <- Call.char $ transposeFromOrder trans
      diagPtr <- Call.char $ charFromTriDiag diag
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.number one
      aPtr <- unpackToTemp (unpack realOrderA) m0 a
      ldaPtr <- Call.leadingDim m0
      betaPtr <- Call.number zero
      bPtr <- ContT $ withForeignPtr b
      ldbPtr <- Call.leadingDim m
      let runTRMM = do
            copyBlock size bPtr cPtr
            BlasGen.trmm sidePtr uploPtr transPtr diagPtr
               mPtr nPtr alphaPtr aPtr ldaPtr cPtr ldbPtr
      liftIO $
         MatrixShape.caseUpLoSym uploA
            runTRMM
            runTRMM
            (BlasGen.symm sidePtr uploPtr
               mPtr nPtr alphaPtr aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldbPtr)


autoUploOrder ::
   (MatrixShape.Content lo, MatrixShape.Content up) => Order -> (Order, (lo,up))
autoUploOrder order =
   case MatrixShape.autoUplo of
      uplo -> (uploOrder uplo order, uplo)
