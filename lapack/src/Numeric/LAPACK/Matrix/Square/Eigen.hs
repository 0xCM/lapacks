{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.Square.Eigen (
   values,
   schur,
   decompose,
   ComplexOf,
   ) where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import Numeric.LAPACK.Matrix.Shape.Private (Order(ColumnMajor), swapOnRowMajor)
import Numeric.LAPACK.Matrix.Private (Square, argSquare)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (ComplexOf, RealOf, zero)
import Numeric.LAPACK.Private
         (copyConjugate, copyToTemp, copyToColumnMajor,
          withAutoWorkspaceInfo)

import qualified Numeric.LAPACK.FFI.Complex as LapackComplex
import qualified Numeric.LAPACK.FFI.Real as LapackReal
import qualified Numeric.BLAS.FFI.Real as BlasReal
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable.Unchecked.Monadic as ArrayIO
import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)

import Foreign.Marshal.Array (advancePtr, peekArray)
import Foreign.C.Types (CInt, CChar)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Ptr (Ptr, nullPtr, nullFunPtr, castPtr)
import Foreign.Storable (Storable)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)

import Data.Complex (Complex)


values ::
   (Shape.C sh, Class.Floating a) =>
   Square sh a -> Vector sh (ComplexOf a)
values =
   getValues $
   Class.switchFloating
      (Values valuesAux) (Values valuesAux)
      (Values valuesAux) (Values valuesAux)

type Values_ sh a = Square sh a -> Vector sh (ComplexOf a)

newtype Values sh a = Values {getValues :: Values_ sh a}

valuesAux ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Storable ar) =>
   Values_ sh a
valuesAux = argSquare $ \_order size a ->
      Array.unsafeCreateWithSize size $ \n wPtr -> do
   let lda = n
   evalContT $ do
      jobvsPtr <- Call.char 'N'
      sortPtr <- Call.char 'N'
      aPtr <- copyToTemp (n*n) a
      ldaPtr <- Call.leadingDim lda
      sdimPtr <- Call.alloca
      let vsPtr = nullPtr
      ldvsPtr <- Call.leadingDim n
      let bworkPtr = nullPtr
      liftIO $
         withAutoWorkspaceInfo eigenMsg "gees" $ \workPtr lworkPtr infoPtr ->
         gees
            jobvsPtr sortPtr n aPtr ldaPtr sdimPtr
            wPtr vsPtr ldvsPtr workPtr lworkPtr bworkPtr infoPtr


{- |
If @(q,r) = schur a@, then @a = q \<#\> r \<#\> adjoint q@,
where @q@ is unitary (orthogonal)
and @r@ is a right-upper triangular matrix for complex @a@
and a 1x1-or-2x2-block upper triangular matrix for real @a@.
With @takeDiagonal r@ you get all eigenvalues of @a@ if @a@ is complex
and the real parts of the eigenvalues if @a@ is real.
Complex conjugated eigenvalues of a real matrix @a@
are encoded as 2x2 blocks along the diagonal.
-}
schur ::
   (Shape.C sh, Class.Floating a) =>
   Square sh a -> (Square sh a, Square sh a)
schur =
   getSchur $
   Class.switchFloating
      (Schur schurAux) (Schur schurAux)
      (Schur schurAux) (Schur schurAux)

type Schur_ sh a = Square sh a -> (Square sh a, Square sh a)

newtype Schur sh a = Schur {getSchur :: Schur_ sh a}

schurAux ::
   (Shape.C sh, Class.Floating a, RealOf a ~ ar, Storable ar) =>
   Schur_ sh a
schurAux = argSquare $ \order size a ->
   let sh = MatrixShape.square ColumnMajor size
   in Array.unsafeCreateWithSizeAndResult sh $ \_ vsPtr ->
      ArrayIO.unsafeCreate sh $ \sPtr -> do

   let n = Shape.size size
   let lda = n
   evalContT $ do
      jobvsPtr <- Call.char 'V'
      sortPtr <- Call.char 'N'
      aPtr <- ContT $ withForeignPtr a
      liftIO $ copyToColumnMajor order n n aPtr sPtr
      ldaPtr <- Call.leadingDim lda
      sdimPtr <- Call.alloca
      wPtr <- Call.allocaArray n
      ldvsPtr <- Call.leadingDim n
      let bworkPtr = nullPtr
      liftIO $
         withAutoWorkspaceInfo eigenMsg "gees" $ \workPtr lworkPtr infoPtr ->
         gees
            jobvsPtr sortPtr n sPtr ldaPtr sdimPtr
            wPtr vsPtr ldvsPtr workPtr lworkPtr bworkPtr infoPtr



type GEES_ ar a =
   Ptr CChar -> Ptr CChar -> Int -> Ptr a -> Ptr CInt ->
   Ptr CInt -> Ptr (Complex ar) -> Ptr a -> Ptr CInt ->
   Ptr a -> Ptr CInt -> Ptr Bool -> Ptr CInt -> IO ()

newtype GEES a = GEES {getGEES :: GEES_ (RealOf a) a}

gees :: Class.Floating a => GEES_ (RealOf a) a
gees =
   getGEES $
   Class.switchFloating
      (GEES geesReal) (GEES geesReal) (GEES geesComplex) (GEES geesComplex)

geesReal :: Class.Real a => GEES_ a a
geesReal
      jobvsPtr sortPtr n aPtr ldaPtr sdimPtr
      wPtr vsPtr ldvsPtr workPtr lworkPtr bworkPtr infoPtr =
   evalContT $ do
      let selectPtr = nullFunPtr
      nPtr <- Call.cint n
      wrPtr <- Call.allocaArray n
      wiPtr <- Call.allocaArray n
      liftIO $
         LapackReal.gees
            jobvsPtr sortPtr selectPtr nPtr aPtr ldaPtr sdimPtr
            wrPtr wiPtr vsPtr ldvsPtr workPtr lworkPtr bworkPtr infoPtr
      liftIO $ zipComplex n wrPtr wiPtr wPtr

geesComplex :: Class.Real a => GEES_ a (Complex a)
geesComplex
      jobvsPtr sortPtr n aPtr ldaPtr sdimPtr
      wPtr vsPtr ldvsPtr workPtr lworkPtr bworkPtr infoPtr =
   evalContT $ do
      let selectPtr = nullFunPtr
      nPtr <- Call.cint n
      rworkPtr <- Call.allocaArray n
      liftIO $
         LapackComplex.gees
            jobvsPtr sortPtr selectPtr nPtr aPtr ldaPtr sdimPtr
            wPtr vsPtr ldvsPtr workPtr lworkPtr rworkPtr bworkPtr infoPtr



decompose ::
   (Shape.C sh, Class.Floating a) =>
   Square sh a ->
   (Square sh (ComplexOf a),
    Vector sh (ComplexOf a),
    Square sh (ComplexOf a))
decompose =
   getDecompose $
   Class.switchFloating
      (Decompose decomposeReal)
      (Decompose decomposeReal)
      (Decompose decomposeComplex)
      (Decompose decomposeComplex)

newtype Decompose sh a =
   Decompose {
      getDecompose ::
         Square sh a ->
         (Square sh (ComplexOf a),
          Vector sh (ComplexOf a),
          Square sh (ComplexOf a))
   }

decomposeReal ::
   (Shape.C sh, Class.Real a) =>
   Square sh a ->
   (Square sh (Complex a), Vector sh (Complex a), Square sh (Complex a))
decomposeReal = argSquare $ \order size a ->
   (\(w, (vlc,vrc)) -> (vlc, w, vrc)) $
   Array.unsafeCreateWithSizeAndResult size $ \n wPtr ->
   evalContT $ do
      jobvlPtr <- Call.char 'V'
      jobvrPtr <- Call.char 'V'
      nPtr <- Call.cint n
      aPtr <- copyToTemp (n*n) a
      ldaPtr <- Call.leadingDim n
      wrPtr <- Call.allocaArray n
      wiPtr <- Call.allocaArray n
      vlPtr <- Call.allocaArray (n*n)
      ldvlPtr <- Call.leadingDim n
      vrPtr <- Call.allocaArray (n*n)
      ldvrPtr <- Call.leadingDim n
      liftIO $ withAutoWorkspaceInfo eigenMsg "geev" $
         LapackReal.geev
            jobvlPtr jobvrPtr nPtr aPtr ldaPtr
            wrPtr wiPtr vlPtr ldvlPtr vrPtr ldvrPtr
      liftIO $ zipComplex n wrPtr wiPtr wPtr
      liftIO $ createArrayPair order (MatrixShape.square ColumnMajor size) $
         \vlcPtr vrcPtr -> do
            eigenvectorsToComplex n wiPtr vlPtr vlcPtr
            eigenvectorsToComplex n wiPtr vrPtr vrcPtr

eigenvectorsToComplex ::
   (Eq a, Class.Real a) =>
   Int -> Ptr a -> Ptr a -> Ptr (Complex a) -> IO ()
eigenvectorsToComplex n wiPtr vPtr vcPtr = evalContT $ do
   nPtr <- Call.cint n
   zeroPtr <- Call.real zero
   inc0Ptr <- Call.cint 0
   inc1Ptr <- Call.cint 1
   inc2Ptr <- Call.cint 2
   liftIO $ do
      let go _ _ [] = return ()
          go xPtr yPtr (False:wi) = do
            let yrPtr = castPtr yPtr
            let yiPtr = advancePtr yrPtr 1
            BlasReal.copy nPtr xPtr    inc1Ptr yrPtr inc2Ptr
            BlasReal.copy nPtr zeroPtr inc0Ptr yiPtr inc2Ptr
            go (advancePtr xPtr n) (advancePtr yPtr n) wi
          go xPtr yPtr (True:True:wi) = do
            let xrPtr = xPtr
            let xiPtr = advancePtr xPtr n
            let yrPtr = castPtr yPtr
            let yiPtr = advancePtr yrPtr 1
            let y1Ptr = advancePtr yPtr n
            BlasReal.copy nPtr xrPtr inc1Ptr yrPtr inc2Ptr
            BlasReal.copy nPtr xiPtr inc1Ptr yiPtr inc2Ptr
            copyConjugate nPtr yPtr inc1Ptr y1Ptr inc1Ptr
            go (advancePtr xPtr (2*n)) (advancePtr yPtr (2*n)) wi
          go _xPtr _yPtr wi =
            error $ "eigenvectorToComplex: invalid non-real pattern " ++ show wi
      go vPtr vcPtr . map (zero/=) =<< peekArray n wiPtr

decomposeComplex ::
   (Shape.C sh, Class.Real a) =>
   Square sh (Complex a) ->
   (Square sh (Complex a), Vector sh (Complex a), Square sh (Complex a))
decomposeComplex = argSquare $ \order size a ->
   (\(w, (vlc,vrc)) -> (vlc, w, vrc)) $
   Array.unsafeCreateWithSizeAndResult size $ \n wPtr ->
   evalContT $ do
      jobvlPtr <- Call.char 'V'
      jobvrPtr <- Call.char 'V'
      nPtr <- Call.cint n
      aPtr <- copyToTemp (n*n) a
      ldaPtr <- Call.leadingDim n
      ldvlPtr <- Call.leadingDim n
      ldvrPtr <- Call.leadingDim n
      rworkPtr <- Call.allocaArray (2*n)

      liftIO $ createArrayPair order (MatrixShape.square ColumnMajor size) $
         \vlPtr vrPtr ->

         withAutoWorkspaceInfo eigenMsg "geev" $ \workPtr lworkPtr infoPtr ->
         LapackComplex.geev
            jobvlPtr jobvrPtr nPtr aPtr ldaPtr
            wPtr vlPtr ldvlPtr vrPtr ldvrPtr
            workPtr lworkPtr rworkPtr infoPtr


zipComplex ::
   (Class.Real a) => Int -> Ptr a -> Ptr a -> Ptr (Complex a) -> IO ()
zipComplex n vr vi vc =
   evalContT $ do
      nPtr <- Call.cint n
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 2
      let yPtr = castPtr vc
      liftIO $ BlasReal.copy nPtr vr incxPtr yPtr incyPtr
      liftIO $ BlasReal.copy nPtr vi incxPtr (advancePtr yPtr 1) incyPtr

createArrayPair ::
   (Shape.C sh, Storable a) =>
   Order -> sh -> (Ptr a -> Ptr a -> IO ()) ->
   IO (Array sh a, Array sh a)
createArrayPair order sh act =
   fmap (swapOnRowMajor order) $
   ArrayIO.unsafeCreateWithSizeAndResult sh $ \_ vrcPtr ->
   ArrayIO.unsafeCreate sh $ \vlcPtr -> act vlcPtr vrcPtr


eigenMsg :: String
eigenMsg = "only eigenvalues starting with the %d-th one converged"
