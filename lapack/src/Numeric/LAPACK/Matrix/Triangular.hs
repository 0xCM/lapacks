{-# LANGUAGE ConstraintKinds #-}
module Numeric.LAPACK.Matrix.Triangular (
   module Numeric.LAPACK.Matrix.Triangular.Basic,
   module Numeric.LAPACK.Matrix.Triangular.Linear,
   size,

   eigenvalues,
   eigensystem,
   ) where

import qualified Numeric.LAPACK.Matrix.Triangular.Eigen as Eigen
import Numeric.LAPACK.Matrix.Triangular.Basic
import Numeric.LAPACK.Matrix.Triangular.Linear

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import Numeric.LAPACK.Matrix.Shape.Private (NonUnit)
import Numeric.LAPACK.Vector (Vector)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape


size :: Triangular lo diag up sh a -> sh
size = MatrixShape.triangularSize . Array.shape


eigenvalues ::
   (MatrixShape.DiagUpLo lo up, Shape.C sh, Class.Floating a) =>
   Triangular lo diag up sh a -> Vector sh a
eigenvalues = Eigen.values


{- |
@(vr,d,vlAdj) = eigensystem a@

Counterintuitively, @vr@ contains the right eigenvectors as columns
and @vlAdj@ contains the left conjugated eigenvectors as rows.
The idea is to provide a decomposition of @a@.
If @a@ is diagonalizable, then @vr@ and @vlAdj@
are almost inverse to each other.
More precisely, @vlAdj \<#\> vr@ is a diagonal matrix.
This is because the eigenvectors are not normalized.
With the following scaling, the decomposition becomes perfect:

> let scal = Array.map recip $ takeDiagonal $ vlAdj <#> vr
> a == vr <#> diagonal d <#> diagonal scal <#> vlAdj

If @a@ is non-diagonalizable
then some columns of @vr@ and corresponding rows of @vlAdj@ are left zero
and the above property does not hold.
-}
eigensystem ::
   (MatrixShape.DiagUpLo lo up, Shape.C sh, Class.Floating a) =>
   Triangular lo NonUnit up sh a ->
   (Triangular lo NonUnit up sh a, Vector sh a, Triangular lo NonUnit up sh a)
eigensystem = Eigen.decompose
