{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.Hermitian (
   module Numeric.LAPACK.Matrix.Hermitian.Basic,
   module Numeric.LAPACK.Matrix.Hermitian.Linear,

   eigenvalues,
   eigensystem,
   ) where

import qualified Numeric.LAPACK.Matrix.Hermitian.Eigen as Eigen
import Numeric.LAPACK.Matrix.Hermitian.Basic
import Numeric.LAPACK.Matrix.Hermitian.Linear

import Numeric.LAPACK.Matrix.Private (Square)
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Shape as Shape


eigenvalues ::
   (Shape.C sh, Class.Floating a) =>
   Hermitian sh a -> Vector sh (RealOf a)
eigenvalues = Eigen.values

{- |
For symmetric eigenvalue problems, @eigensystem@ and @schur@ coincide.
-}
eigensystem ::
   (Shape.C sh, Class.Floating a) =>
   Hermitian sh a -> (Square sh a, Vector sh (RealOf a))
eigensystem = Eigen.decompose
