{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.Hermitian.Private where

import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf)


newtype Diagonal f sh a =
   Diagonal {runDiagonal :: Vector sh (RealOf a) -> f a}

newtype TakeDiagonal f sh a =
   TakeDiagonal {runTakeDiagonal :: f a -> Vector sh (RealOf a)}

newtype Determinant f a = Determinant {getDeterminant :: f a -> RealOf a}
