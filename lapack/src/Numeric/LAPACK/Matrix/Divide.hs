{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE UndecidableInstances #-}
module Numeric.LAPACK.Matrix.Divide where

import qualified Numeric.LAPACK.Matrix.Square.Linear
                                           as Square
import qualified Numeric.LAPACK.Matrix.Triangular.Linear
                                           as Triangular
import qualified Numeric.LAPACK.Matrix.Hermitian.Linear
                                           as Hermitian
import qualified Numeric.LAPACK.Matrix.Banded.Linear
                                           as Banded
import qualified Numeric.LAPACK.Matrix.BandedHermitianPositiveDefinite.Linear
                                           as BandedHermitianPositiveDefinite

import qualified Numeric.LAPACK.Matrix.Basic as Basic
import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Shape.Private (HeightOf)
import Numeric.LAPACK.Matrix.Extent.Private (Small)
import Numeric.LAPACK.Matrix.Private (Full)
import Numeric.LAPACK.Vector (Vector)

import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary as Unary

import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)


class (Shape.C shape) => Solve shape where
   solve ::
      (Class.Floating a, HeightOf shape ~ height, Eq height,
       Extent.C horiz, Extent.C vert, Shape.C nrhs) =>
      Array shape a ->
      Full vert horiz height nrhs a -> Full vert horiz height nrhs a

class (Solve shape) => Inverse shape where
   inverse :: (Class.Floating a) => Array shape a -> Array shape a

solveVector ::
   (Solve shape, HeightOf shape ~ height, Eq height, Class.Floating a) =>
   Array shape a -> Vector height a -> Vector height a
solveVector m =
   Basic.flattenColumn . solve m . Basic.singleColumn MatrixShape.ColumnMajor


instance
   (vert ~ Small, horiz ~ Small,
    Shape.C width, Shape.C height, height ~ width) =>
      Solve (MatrixShape.Full vert horiz height width) where
   solve = Square.solve

instance
   (vert ~ Small, horiz ~ Small,
    Shape.C width, Shape.C height, height ~ width) =>
      Inverse (MatrixShape.Full vert horiz height width) where
   inverse = Square.inverse


instance (Shape.C shape) => Solve (MatrixShape.Hermitian shape) where
   solve = Hermitian.solve

instance (Shape.C shape) => Inverse (MatrixShape.Hermitian shape) where
   inverse = Hermitian.inverse


instance
   (MatrixShape.Content lo, MatrixShape.Content up,
    MatrixShape.TriDiag diag, Shape.C shape) =>
      Solve (MatrixShape.Triangular lo diag up shape) where
   solve = Triangular.solve

instance
   (MatrixShape.DiagUpLo lo up,
    MatrixShape.TriDiag diag, Shape.C shape) =>
      Inverse (MatrixShape.Triangular lo diag up shape) where
   inverse = Triangular.inverse


instance
   (Unary.Natural sub, Unary.Natural super, vert ~ Small, horiz ~ Small,
    Shape.C width, Shape.C height, width ~ height) =>
      Solve (MatrixShape.Banded sub super vert horiz height width) where
   solve = Banded.solve


{- |
There is no solver for indefinite matrices.
Thus the instance will fail for indefinite but solvable systems.
-}
instance
   (Unary.Natural offDiag, Shape.C size) =>
      Solve (MatrixShape.BandedHermitian offDiag size) where
   solve = BandedHermitianPositiveDefinite.solve
