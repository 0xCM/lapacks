module Numeric.LAPACK.Matrix.Square (
   module Numeric.LAPACK.Matrix.Square.Basic,
   module Numeric.LAPACK.Matrix.Square.Linear,

   eigenvalues,
   Eigen.schur,
   eigensystem,
   ComplexOf,
   ) where

import qualified Numeric.LAPACK.Matrix.Square.Eigen as Eigen
import Numeric.LAPACK.Matrix.Square.Basic
import Numeric.LAPACK.Matrix.Square.Linear

import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (ComplexOf)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Shape as Shape


eigenvalues ::
   (Shape.C sh, Class.Floating a) =>
   Square sh a -> Vector sh (ComplexOf a)
eigenvalues = Eigen.values

{- |
@(vr,d,vl) = eigensystem a@

Counterintuitively, @vr@ contains the right eigenvectors
and @vl@ contains the left eigenvectors as columns.
The idea is to provide a decomposition of @a@.
If @a@ is diagonalizable, then @vr@ and @vl@ are almost inverse to each other.
More precisely, @adjoint vl \<#\> vr@ is a diagonal matrix.
This is because all eigenvectors are normalized to Euclidean norm 1.
With the following scaling, the decomposition becomes perfect:

> let scal = Array.map recip $ takeDiagonal $ adjoint vl <#> vr
> a == vr <#> diagonal d <#> diagonal scal <#> adjoint vl

If @a@ is non-diagonalizable then some columns of @vr@ and @vl@ are left zero
and the above property does not hold.
-}
eigensystem ::
   (Shape.C sh, Class.Floating a) =>
   Square sh a ->
   (Square sh (ComplexOf a),
    Vector sh (ComplexOf a),
    Square sh (ComplexOf a))
eigensystem = Eigen.decompose
