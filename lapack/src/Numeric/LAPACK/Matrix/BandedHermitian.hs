{-# LANGUAGE TypeFamilies #-}
module Numeric.LAPACK.Matrix.BandedHermitian (
   module Numeric.LAPACK.Matrix.BandedHermitian.Basic,

   eigenvalues,
   eigensystem,
   ) where

import qualified Numeric.LAPACK.Matrix.BandedHermitian.Eigen as Eigen
import Numeric.LAPACK.Matrix.BandedHermitian.Basic

import qualified Numeric.LAPACK.Matrix.Private as Matrix
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf)

import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary as Unary

import qualified Data.Array.Comfort.Shape as Shape


eigenvalues ::
   (Unary.Natural offDiag, Shape.C sh, Class.Floating a) =>
   BandedHermitian offDiag sh a -> Vector sh (RealOf a)
eigenvalues = Eigen.values

{- |
For symmetric eigenvalue problems, @eigensystem@ and @schur@ coincide.
-}
eigensystem ::
   (Unary.Natural offDiag, Shape.C sh, Class.Floating a) =>
   BandedHermitian offDiag sh a -> (Matrix.Square sh a, Vector sh (RealOf a))
eigensystem = Eigen.decompose
