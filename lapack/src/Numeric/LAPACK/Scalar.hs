{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Scalar (
   RealOf,
   ComplexOf,
   zero,
   one,
   minusOne,
   isZero,
   selectReal,
   selectFloating,

   fromReal,
   absolute,
   absoluteSquared,
   norm1,
   realPart,
   conjugate,
   ) where

import Numeric.LAPACK.Wrapper (Flip(Flip, getFlip))

import qualified Numeric.Netlib.Class as Class

import Data.Functor.Identity (Identity(Identity, runIdentity))

import qualified Data.Complex as Complex
import Data.Complex (Complex((:+)))
import Data.Monoid (Endo(Endo,appEndo))


type family RealOf x

type instance RealOf Float = Float
type instance RealOf Double = Double
type instance RealOf (Complex a) = a


type ComplexOf x = Complex (RealOf x)


-- move to netlib-carray:Utility or netlib-ffi:Class
zero, one, minusOne :: Class.Floating a => a
zero = selectFloating 0 0 0 0
one = selectFloating 1 1 1 1
minusOne = selectFloating (-1) (-1) (-1) (-1)

selectReal :: (Class.Real a) => Float -> Double -> a
selectReal rf rd =
   runIdentity $ Class.switchReal (Identity rf) (Identity rd)

selectFloating ::
   (Class.Floating a) =>
   Float -> Double -> Complex Float -> Complex Double -> a
selectFloating rf rd cf cd =
   runIdentity $
   Class.switchFloating
      (Identity rf) (Identity rd) (Identity cf) (Identity cd)


isZero :: Class.Floating a => a -> Bool
isZero =
   getFlip $
   Class.switchFloating
      (Flip (0==)) (Flip (0==))
      (Flip (0==)) (Flip (0==))


newtype FromReal a = FromReal {getFromReal :: RealOf a -> a}

fromReal :: (Class.Floating a) => RealOf a -> a
fromReal =
   getFromReal $
   Class.switchFloating
      (FromReal id)
      (FromReal id)
      (FromReal (:+0))
      (FromReal (:+0))

newtype ToReal a = ToReal {getToReal :: a -> RealOf a}

realPart :: (Class.Floating a) => a -> RealOf a
realPart =
   getToReal $
   Class.switchFloating
      (ToReal id)
      (ToReal id)
      (ToReal Complex.realPart)
      (ToReal Complex.realPart)

absolute :: (Class.Floating a) => a -> RealOf a
absolute =
   getToReal $
   Class.switchFloating
      (ToReal abs)
      (ToReal abs)
      (ToReal Complex.magnitude)
      (ToReal Complex.magnitude)


norm1 :: (Class.Floating a) => a -> RealOf a
norm1 =
   getToReal $
   Class.switchFloating
      (ToReal abs)
      (ToReal abs)
      (ToReal norm1Complex)
      (ToReal norm1Complex)

norm1Complex :: (Class.Real a) => Complex a -> a
norm1Complex (r:+i) = abs r + abs i


absoluteSquared :: (Class.Floating a) => a -> RealOf a
absoluteSquared =
   getToReal $
   Class.switchFloating
      (ToReal absoluteSquaredReal)
      (ToReal absoluteSquaredReal)
      (ToReal absoluteSquaredComplex)
      (ToReal absoluteSquaredComplex)

absoluteSquaredReal :: (Class.Real a) => a -> a
absoluteSquaredReal a = a*a

absoluteSquaredComplex :: (Class.Real a) => Complex a -> a
absoluteSquaredComplex (r:+i) = r*r+i*i


conjugate :: (Class.Floating a) => a -> a
conjugate =
   appEndo $
   Class.switchFloating
      (Endo id) (Endo id) (Endo Complex.conjugate) (Endo Complex.conjugate)
