{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Numeric.LAPACK.Vector (
   Vector,
   RealOf,
   ComplexOf,
   toList,
   fromList,
   autoFromList,
   append, take, drop,
   takeLeft, takeRight,
   constant,
   unit,
   dot, inner,
   sum,
   absSum,
   norm1,
   norm2,
   normInf,
   normInf1,
   argAbsMaximum,
   argAbs1Maximum,
   product,
   scale, scaleReal,
   add, sub,
   mac,
   mul,

   conjugate,
   fromReal,
   toComplex,
   realPart,
   complexFromReal,
   complexToRealPart,
   complexToImaginaryPart,
   zipComplex,
   unzipComplex,

   random, RandomDistribution(..),
   ) where

import qualified Numeric.LAPACK.Scalar as Scalar
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Matrix.Private (ZeroInt)
import Numeric.LAPACK.Scalar (ComplexOf, RealOf, zero, one, minusOne, absolute)
import Numeric.LAPACK.Private (fill, copyConjugate)

import qualified Numeric.LAPACK.FFI.Generic as LapackGen
import qualified Numeric.LAPACK.FFI.Complex as LapackComplex
import qualified Numeric.BLAS.FFI.Generic as BlasGen
import qualified Numeric.BLAS.FFI.Complex as BlasComplex
import qualified Numeric.BLAS.FFI.Real as BlasReal
import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import Foreign.Marshal.Array (copyArray, advancePtr)
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Ptr (Ptr, castPtr)
import Foreign.Storable (Storable, peek, peekElemOff, pokeElemOff)
import Foreign.C.Types (CInt)

import System.IO.Unsafe (unsafePerformIO)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)
import Control.Applicative (Const(Const,getConst), liftA3, (<$>))

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))
import Data.Array.Comfort.Shape ((:+:)((:+:)))

import Data.Complex (Complex)
import Data.Tuple.HT (mapFst, uncurry3)
import Data.Word (Word64)
import Data.Bits (shiftR, (.&.))

import Prelude hiding (sum, product, take, drop)


type Vector = Array


toList :: (Shape.C sh, Storable a) => Vector sh a -> [a]
toList = Array.toList

fromList :: (Shape.C sh, Storable a) => sh -> [a] -> Vector sh a
fromList = Array.fromList

autoFromList :: (Storable a) => [a] -> Vector (Shape.ZeroBased Int) a
autoFromList = Array.vectorFromList


constant :: (Shape.C sh, Class.Floating a) => sh -> a -> Vector sh a
constant sh a = Array.unsafeCreateWithSize sh $ fill a

unit ::
   (Shape.Indexed sh, Class.Floating a) =>
   sh -> Shape.Index sh -> Vector sh a
unit sh ix = Array.unsafeCreateWithSize sh $ \n xPtr -> do
   fill zero n xPtr
   pokeElemOff xPtr (Shape.offset sh ix) one


append ::
   (Shape.C shx, Shape.C shy, Storable a) =>
   Vector shx a -> Vector shy a -> Vector (shx:+:shy) a
append (Array shX x) (Array shY y) =
   Array.unsafeCreate (shX:+:shY) $ \zPtr ->
   evalContT $ do
      xPtr <- ContT $ withForeignPtr x
      yPtr <- ContT $ withForeignPtr y
      let sizeX = Shape.size shX
      let sizeY = Shape.size shY
      liftIO $ do
         copyArray zPtr xPtr sizeX
         copyArray (advancePtr zPtr sizeX) yPtr sizeY

take, drop :: (Storable a) => Int -> Vector ZeroInt a -> Vector ZeroInt a
take n = takeLeft . split n
drop n = takeRight . split n

split :: (Storable a) => Int -> Vector ZeroInt a -> Vector (ZeroInt:+:ZeroInt) a
split n =
   Array.mapShape
      (\(Shape.ZeroBased m) ->
         if n<0
            then error "Vector.split: negative number of elements"
            else
               let k = min n m
               in Shape.ZeroBased k :+: Shape.ZeroBased (m-k))

takeLeft ::
   (Shape.C sh0, Shape.C sh1, Storable a) =>
   Vector (sh0:+:sh1) a -> Vector sh0 a
takeLeft (Array (sh0 :+: _sh1) x) =
   Array.unsafeCreateWithSize sh0 $ \k yPtr ->
   withForeignPtr x $ \xPtr -> copyArray yPtr xPtr k

takeRight ::
   (Shape.C sh0, Shape.C sh1, Storable a) =>
   Vector (sh0:+:sh1) a -> Vector sh1 a
takeRight (Array (sh0:+:sh1) x) =
   Array.unsafeCreateWithSize sh1 $ \k yPtr ->
   withForeignPtr x $ \xPtr ->
      copyArray yPtr (advancePtr xPtr (Shape.size sh0)) k


newtype Dot sh a = Dot {runDot :: Vector sh a -> Vector sh a -> a}

{- |
> dot x y = Matrix.toScalar (singleRow x <#> singleColumn y)
-}
dot ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Vector sh a -> Vector sh a -> a
dot =
   runDot $
   Class.switchFloating
      (Dot dotReal)
      (Dot dotReal)
      (Dot $ dotComplex 'T')
      (Dot $ dotComplex 'T')

{- |
> inner x y = dot (conjugate x) y
-}
inner ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Vector sh a -> Vector sh a -> a
inner =
   runDot $
   Class.switchFloating
      (Dot dotReal)
      (Dot dotReal)
      (Dot $ dotComplex 'C')
      (Dot $ dotComplex 'C')

dotReal ::
   (Shape.C sh, Eq sh, Class.Real a) =>
   Vector sh a -> Vector sh a -> a
dotReal arrX@(Array shX _x) (Array shY y) = unsafePerformIO $ do
   Call.assert "dot: shapes mismatch" (shX == shY)
   evalContT $ do
      (nPtr, sxPtr, incxPtr) <- vectorArgs arrX
      syPtr <- ContT $ withForeignPtr y
      incyPtr <- Call.cint 1
      liftIO $ BlasReal.dot nPtr sxPtr incxPtr syPtr incyPtr

{-
We cannot use 'cdot' because Haskell's FFI
does not support Complex numbers as return values.
-}
dotComplex ::
   (Shape.C sh, Eq sh, Class.Real a) =>
   Char -> Vector sh (Complex a) -> Vector sh (Complex a) -> Complex a
dotComplex trans (Array shX x) (Array shY y) = unsafePerformIO $ do
   Call.assert "dot: shapes mismatch" (shX == shY)
   evalContT $ do
      let m = Shape.size shX
      transPtr <- Call.char trans
      mPtr <- Call.cint m
      nPtr <- Call.cint 1
      alphaPtr <- Call.number one
      xPtr <- ContT $ withForeignPtr x
      ldxPtr <- Call.leadingDim m
      yPtr <- ContT $ withForeignPtr y
      incyPtr <- Call.cint 1
      betaPtr <- Call.number zero
      zPtr <- Call.alloca
      inczPtr <- Call.cint 1
      liftIO $
         Private.gemv
            transPtr mPtr nPtr alphaPtr xPtr ldxPtr
            yPtr incyPtr betaPtr zPtr inczPtr
      liftIO $ peek zPtr

sum :: (Shape.C sh, Class.Floating a) => Vector sh a -> a
sum (Array sh x) = unsafePerformIO $
   withForeignPtr x $ \xPtr -> Private.sum (Shape.size sh) xPtr 1

norm1 :: (Shape.C sh, Class.Floating a) => Vector sh a -> RealOf a
norm1 arr = unsafePerformIO $
   evalContT $ liftIO . uncurry3 csum1 =<< vectorArgs arr

csum1 :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr CInt -> IO (RealOf a)
csum1 =
   getNorm $
   Class.switchFloating
      (Norm BlasReal.asum)
      (Norm BlasReal.asum)
      (Norm LapackComplex.sum1)
      (Norm LapackComplex.sum1)


{- |
Sum of the absolute values of real numbers or components of complex numbers.
For real numbers it is equivalent to 'norm1'.
-}
absSum :: (Shape.C sh, Class.Floating a) => Vector sh a -> RealOf a
absSum arr = unsafePerformIO $
   evalContT $ liftIO . uncurry3 asum =<< vectorArgs arr

asum :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr CInt -> IO (RealOf a)
asum =
   getNorm $
   Class.switchFloating
      (Norm BlasReal.asum) (Norm BlasReal.asum)
      (Norm BlasComplex.casum) (Norm BlasComplex.casum)


{- |
Euclidean norm of a vector or Frobenius norm of a matrix.
-}
norm2 :: (Shape.C sh, Class.Floating a) => Vector sh a -> RealOf a
norm2 arr = unsafePerformIO $
   evalContT $ liftIO . uncurry3 nrm2 =<< vectorArgs arr

nrm2 :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr CInt -> IO (RealOf a)
nrm2 =
   getNorm $
   Class.switchFloating
      (Norm BlasReal.nrm2) (Norm BlasReal.nrm2)
      (Norm BlasComplex.cnrm2) (Norm BlasComplex.cnrm2)

newtype Norm a =
   Norm {getNorm :: Ptr CInt -> Ptr a -> Ptr CInt -> IO (RealOf a)}


normInf :: (Shape.C sh, Class.Floating a) => Vector sh a -> RealOf a
normInf arr = unsafePerformIO $
   evalContT $ do
      (nPtr, sxPtr, incxPtr) <- vectorArgs arr
      liftIO $
         fmap (absolute . maybe zero snd) $
         peekElemOff1 sxPtr =<< absMax nPtr sxPtr incxPtr

{- |
Computes (almost) the infinity norm of the vector.
For complex numbers every element is replaced
by the sum of the absolute component values first.
-}
normInf1 :: (Shape.C sh, Class.Floating a) => Vector sh a -> RealOf a
normInf1 arr = unsafePerformIO $
   evalContT $ do
      (nPtr, sxPtr, incxPtr) <- vectorArgs arr
      liftIO $
         fmap (Scalar.norm1 . maybe zero snd) $
         peekElemOff1 sxPtr =<< BlasGen.iamax nPtr sxPtr incxPtr


{- |
Returns the index and value of the element with the maximal absolute value.
Caution: It actually returns the value of the element, not its absolute value!
-}
argAbsMaximum ::
   (Shape.InvIndexed sh, Class.Floating a) =>
   Vector sh a -> (Shape.Index sh, a)
argAbsMaximum arr = unsafePerformIO $
   evalContT $ do
      (nPtr, sxPtr, incxPtr) <- vectorArgs arr
      liftIO $
         fmap
            (maybe
               (error "Vector.argAbsMaximum: empty vector")
               (mapFst (Shape.uncheckedIndexFromOffset $ Array.shape arr))) $
         peekElemOff1 sxPtr =<< absMax nPtr sxPtr incxPtr

newtype ArgMaximum a =
   ArgMaximum {runArgMaximum :: Ptr CInt -> Ptr a -> Ptr CInt -> IO CInt}

absMax :: Class.Floating a => Ptr CInt -> Ptr a -> Ptr CInt -> IO CInt
absMax =
   runArgMaximum $
   Class.switchFloating
      (ArgMaximum BlasGen.iamax)
      (ArgMaximum BlasGen.iamax)
      (ArgMaximum LapackComplex.imax1)
      (ArgMaximum LapackComplex.imax1)


{- |
Returns the index and value of the element with the maximal absolute value.
The function does not strictly compare the absolute value of a complex number
but the sum of the absolute complex components.
Caution: It actually returns the value of the element, not its absolute value!
-}
argAbs1Maximum ::
   (Shape.InvIndexed sh, Class.Floating a) =>
   Vector sh a -> (Shape.Index sh, a)
argAbs1Maximum arr = unsafePerformIO $
   evalContT $ do
      (nPtr, sxPtr, incxPtr) <- vectorArgs arr
      liftIO $
         fmap
            (maybe
               (error "Vector.argAbs1Maximum: empty vector")
               (mapFst (Shape.uncheckedIndexFromOffset $ Array.shape arr))) $
         peekElemOff1 sxPtr =<< BlasGen.iamax nPtr sxPtr incxPtr

vectorArgs ::
   (Shape.C sh) => Array sh a -> ContT r IO (Ptr CInt, Ptr a, Ptr CInt)
vectorArgs (Array sh x) =
   liftA3 (,,)
      (Call.cint $ Shape.size sh)
      (ContT $ withForeignPtr x)
      (Call.cint 1)

peekElemOff1 :: (Storable a) => Ptr a -> CInt -> IO (Maybe (Int, a))
peekElemOff1 ptr k1 =
   let k1i = fromIntegral k1
       ki = k1i-1
   in if k1i == 0
         then return Nothing
         else Just . (,) ki <$> peekElemOff ptr ki


product :: (Shape.C sh, Class.Floating a) => Vector sh a -> a
product (Array sh x) = unsafePerformIO $
   withForeignPtr x $ \xPtr -> Private.product (Shape.size sh) xPtr 1


scale, _scale ::
   (Shape.C sh, Class.Floating a) =>
   a -> Vector sh a -> Vector sh a
scale alpha (Array sh x) = Array.unsafeCreateWithSize sh $ \n syPtr -> do
   evalContT $ do
      alphaPtr <- Call.number alpha
      nPtr <- Call.cint n
      sxPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 1
      liftIO $ BlasGen.copy nPtr sxPtr incxPtr syPtr incyPtr
      liftIO $ BlasGen.scal nPtr alphaPtr syPtr incyPtr

_scale a (Array sh b) = Array.unsafeCreateWithSize sh $ \n cPtr -> do
   let m = 1
   let k = 1
   evalContT $ do
      transaPtr <- Call.char 'N'
      transbPtr <- Call.char 'N'
      mPtr <- Call.cint m
      kPtr <- Call.cint k
      nPtr <- Call.cint n
      alphaPtr <- Call.number one
      aPtr <- Call.number a
      ldaPtr <- Call.leadingDim m
      bPtr <- ContT $ withForeignPtr b
      ldbPtr <- Call.leadingDim k
      betaPtr <- Call.number zero
      ldcPtr <- Call.leadingDim m
      liftIO $
         BlasGen.gemm
            transaPtr transbPtr mPtr nPtr kPtr alphaPtr
            aPtr ldaPtr bPtr ldbPtr betaPtr cPtr ldcPtr


scaleReal ::
   (Shape.C sh, Class.Floating a) =>
   RealOf a -> Vector sh a -> Vector sh a
scaleReal =
   getScaleReal $
   Class.switchFloating
      (ScaleReal scale)
      (ScaleReal scale)
      (ScaleReal scaleRealComplex)
      (ScaleReal scaleRealComplex)

newtype ScaleReal f a = ScaleReal {getScaleReal :: RealOf a -> f a -> f a}

scaleRealComplex ::
   (Shape.C sh, Class.Real a) =>
   a -> Vector sh (Complex a) -> Vector sh (Complex a)
scaleRealComplex alpha (Array sh x) =
      Array.unsafeCreateWithSize sh $ \n cyPtr ->
   evalContT $ do
      alphaPtr <- Call.number alpha
      n2Ptr <- Call.cint (2*n)
      cxPtr <- ContT $ withForeignPtr x
      let sxPtr = castPtr cxPtr
      let syPtr = castPtr cyPtr
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 1
      liftIO $ do
         BlasReal.copy n2Ptr sxPtr incxPtr syPtr incyPtr
         BlasReal.scal n2Ptr alphaPtr syPtr incyPtr


add, sub ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Vector sh a -> Vector sh a -> Vector sh a
add = mac one
sub x y = mac minusOne y x

mac ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   a -> Vector sh a -> Vector sh a -> Vector sh a
mac alpha (Array shX x) (Array shY y) =
      Array.unsafeCreateWithSize shX $ \n szPtr -> do
   Call.assert "mac: shapes mismatch" (shX == shY)
   evalContT $ do
      nPtr <- Call.cint n
      saPtr <- Call.number alpha
      sxPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint 1
      syPtr <- ContT $ withForeignPtr y
      incyPtr <- Call.cint 1
      inczPtr <- Call.cint 1
      liftIO $ BlasGen.copy nPtr syPtr incyPtr szPtr inczPtr
      liftIO $ BlasGen.axpy nPtr saPtr sxPtr incxPtr szPtr inczPtr

mul ::
   (Shape.C sh, Eq sh, Class.Floating a) =>
   Vector sh a -> Vector sh a -> Vector sh a
mul (Array shA a) (Array shX x) =
      Array.unsafeCreateWithSize shX $ \n yPtr -> do
   Call.assert "mul: shapes mismatch" (shA == shX)
   evalContT $ do
      transPtr <- Call.char 'N'
      nPtr <- Call.cint n
      klPtr <- Call.cint 0
      kuPtr <- Call.cint 0
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim 1
      xPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint 1
      betaPtr <- Call.number zero
      incyPtr <- Call.cint 1
      liftIO $
         BlasGen.gbmv transPtr
            nPtr nPtr klPtr kuPtr alphaPtr aPtr ldaPtr
            xPtr incxPtr betaPtr yPtr incyPtr


newtype Conjugate sh a = Conjugate {getConjugate :: Vector sh a -> Vector sh a}

conjugate ::
   (Shape.C sh, Class.Floating a) =>
   Vector sh a -> Vector sh a
conjugate =
   getConjugate $
   Class.switchFloating
      (Conjugate id)
      (Conjugate id)
      (Conjugate complexConjugate)
      (Conjugate complexConjugate)

complexConjugate ::
   (Shape.C sh, Class.Real a) =>
   Vector sh (Complex a) -> Vector sh (Complex a)
complexConjugate (Array sh x) = Array.unsafeCreateWithSize sh $ \n syPtr ->
   evalContT $ do
      nPtr <- Call.cint n
      sxPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 1
      liftIO $ copyConjugate nPtr sxPtr incxPtr syPtr incyPtr


fromReal ::
   (Shape.C sh, Class.Floating a) => Vector sh (RealOf a) -> Vector sh a
fromReal =
   getFromReal $
   Class.switchFloating
      (FromReal id)
      (FromReal id)
      (FromReal complexFromReal)
      (FromReal complexFromReal)

newtype FromReal f a = FromReal {getFromReal :: f (RealOf a) -> f a}

toComplex ::
   (Shape.C sh, Class.Floating a) => Vector sh a -> Vector sh (ComplexOf a)
toComplex =
   getToComplex $
   Class.switchFloating
      (ToComplex complexFromReal)
      (ToComplex complexFromReal)
      (ToComplex id)
      (ToComplex id)

newtype ToComplex f a = ToComplex {getToComplex :: f a -> f (ComplexOf a)}

complexFromReal ::
   (Shape.C sh, Class.Real a) => Vector sh a -> Vector sh (Complex a)
complexFromReal (Array sh x) =
   Array.unsafeCreateWithSize sh $ \n yPtr ->
   case castPtr yPtr of
      yrPtr -> evalContT $ do
         nPtr <- Call.cint n
         xPtr <- ContT $ withForeignPtr x
         incxPtr <- Call.cint 1
         incyPtr <- Call.cint 2
         inczPtr <- Call.cint 0
         zPtr <- Call.number zero
         liftIO $ do
            BlasGen.copy nPtr xPtr incxPtr yrPtr incyPtr
            BlasGen.copy nPtr zPtr inczPtr (advancePtr yrPtr 1) incyPtr


realPart ::
   (Shape.C sh, Class.Floating a) => Vector sh a -> Vector sh (RealOf a)
realPart =
   getToReal $
   Class.switchFloating
      (ToReal id)
      (ToReal id)
      (ToReal complexToRealPart)
      (ToReal complexToRealPart)

newtype ToReal f a = ToReal {getToReal :: f a -> f (RealOf a)}


zipComplex ::
   (Shape.C sh, Eq sh, Class.Real a) =>
   Vector sh a -> Vector sh a -> Vector sh (Complex a)
zipComplex (Array shr xr) (Array shi xi) =
   Array.unsafeCreateWithSize shr $ \n yPtr -> evalContT $ do
      liftIO $ Call.assert "zipComplex: shapes mismatch" (shr==shi)
      nPtr <- Call.cint n
      xrPtr <- ContT $ withForeignPtr xr
      xiPtr <- ContT $ withForeignPtr xi
      let yrPtr = castPtr yPtr
      incxPtr <- Call.cint 1
      incyPtr <- Call.cint 2
      liftIO $ do
         BlasGen.copy nPtr xrPtr incxPtr yrPtr incyPtr
         BlasGen.copy nPtr xiPtr incxPtr (advancePtr yrPtr 1) incyPtr


complexToRealPart, complexToImaginaryPart ::
   (Shape.C sh, Class.Real a) => Vector sh (Complex a) -> Vector sh a
complexToRealPart = complexToPart 0
complexToImaginaryPart = complexToPart 1

complexToPart ::
   (Shape.C sh, Class.Real a) => Int -> Vector sh (Complex a) -> Vector sh a
complexToPart offset (Array sh x) =
   Array.unsafeCreateWithSize sh $ \n yPtr -> evalContT $ do
      nPtr <- Call.cint n
      xPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint 2
      incyPtr <- Call.cint 1
      liftIO $
         BlasGen.copy nPtr
            (advancePtr (castPtr xPtr) offset) incxPtr yPtr incyPtr

unzipComplex ::
   (Shape.C sh, Class.Real a) =>
   Vector sh (Complex a) -> (Vector sh a, Vector sh a)
unzipComplex x = (complexToRealPart x, complexToImaginaryPart x)


data RandomDistribution =
     UniformBox01
   | UniformBoxPM1
   | Normal
   | UniformDisc
   | UniformCircle
   deriving (Eq, Ord, Show, Enum)

{-
@random distribution shape seed@

Only the least significant 47 bits of @seed@ are used.
-}
random ::
   (Shape.C sh, Class.Floating a) =>
   RandomDistribution -> sh -> Word64 -> Vector sh a
random dist sh seed = Array.unsafeCreateWithSize sh $ \n xPtr ->
   evalContT $ do
      nPtr <- Call.cint n
      distPtr <-
         Call.cint $
         case (getConst $ isComplexInFunctor xPtr, dist) of
            (_, UniformBox01) -> 1
            (_, UniformBoxPM1) -> 2
            (_, Normal) -> 3
            (True, UniformDisc) -> 4
            (True, UniformCircle) -> 5
            (False, UniformDisc) -> 2
            (False, UniformCircle) ->
               error
                  "Vector.random: UniformCircle not supported for real numbers"
      iseedPtr <- Call.allocaArray 4
      liftIO $ do
         pokeElemOff iseedPtr 0 $ fromIntegral ((seed `shiftR` 35) .&. 0xFFF)
         pokeElemOff iseedPtr 1 $ fromIntegral ((seed `shiftR` 23) .&. 0xFFF)
         pokeElemOff iseedPtr 2 $ fromIntegral ((seed `shiftR` 11) .&. 0xFFF)
         pokeElemOff iseedPtr 3 $ fromIntegral ((seed.&.0x7FF)*2+1)
         LapackGen.larnv distPtr iseedPtr nPtr xPtr

isComplexInFunctor :: (Class.Floating a) => f a -> Const Bool a
isComplexInFunctor _ =
   Class.switchFloating (Const False) (Const False) (Const True) (Const True)
