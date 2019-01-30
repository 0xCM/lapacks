{-# LANGUAGE TypeFamilies #-}
module Test.Square (testsVar) where

import qualified Test.Generator as Gen
import qualified Test.Utility as Util
import Test.Generator ((<|*|>), (<|\|>))
import Test.Utility (approx, approxArray, approxArrayTol, approxMatrix, Tagged)

import qualified Numeric.LAPACK.Matrix.Triangular as Tri
import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix.Square (Square)
import Numeric.LAPACK.Matrix (ZeroInt, (<#>))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, absolute, selectReal)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as Array

import Control.Applicative ((<$>))

import Data.Function.HT (nest)

import qualified Test.QuickCheck as QC


multiplySquare ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
multiplySquare a =
   approxArray (Square.square a) (Square.multiply a a)

multiplyPower ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Int -> Square ZeroInt a -> Bool
multiplyPower n a =
   let b = Square.power (fromIntegral n) a
       c = nest n (Square.multiply a) $ Square.identityFrom a
   in approxArrayTol (1e-6 * (Vector.normInf1 b + Vector.normInf1 c)) b c


determinantSingleton ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   a -> Bool
determinantSingleton a =
   approx 1e-5 a (Square.determinant $ Square.autoFromList [a])

determinantTranspose ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
determinantTranspose a =
   approx 1e-5
      (Square.determinant a) (Square.determinant $ Square.transpose a)


multiplyDeterminant ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Square ZeroInt a, Square ZeroInt a) -> Bool
multiplyDeterminant (a,b) =
   let detA = Square.determinant a
       detB = Square.determinant b
   in approx
         (1e-2 * max 1 (absolute detA) * max 1 (absolute detB))
         (Square.determinant (a<#>b))
         (detA * detB)

multiplyInverse ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
multiplyInverse a =
   let eye = Square.inverse a <#> a
   in approxArrayTol 1e-4 eye (Square.identityFrom eye)


multiplySolve ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Square ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
multiplySolve (a, b) =
   approxMatrix 1e-2 (a <#> Square.solve a b) b

schur ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
schur a =
   let (q,r) = Square.schur a
   in  approxMatrix 1e-4 a (q <#> r <#> Square.adjoint q)


diagonal :: (Class.Floating a) => Vector ZeroInt a -> Tri.Diagonal ZeroInt a
diagonal = Tri.diagonal MatrixShape.ColumnMajor

genDiagonalizable ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Gen.Matrix a Int Int (Square ZeroInt a)
genDiagonalizable = flip Gen.mapGen Gen.invertible $ \ _maxElem a -> do
   d <- Util.genDistinct 3 10 (Square.size a)
   return $ Square.solve a $ diagonal d <#> a


eigensystem ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
eigensystem a =
   let (vr,d,vl) = Square.eigensystem a
       scal = Array.map recip $ Square.takeDiagonal $ Square.adjoint vl <#> vr
   in  approxMatrix (selectReal 1e-1 1e-5)
         (Vector.toComplex a)
         (vr <#> diagonal d <#> diagonal scal <#> Square.adjoint vl)

eigensystemLeft ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
eigensystemLeft a =
   let (_vr,d,vl) = Square.eigensystem a
       vlAdj = Square.adjoint vl
   in  approxMatrix (selectReal 1e-1 1e-5)
         (Vector.toComplex a)
         (Square.solve vlAdj $ diagonal d <#> vlAdj)

eigensystemRight ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Square ZeroInt a -> Bool
eigensystemRight a =
   let (vr,d,_vl) = Square.eigensystem a
       solveLeft b m =
         Matrix.transpose $
         Square.solve (Matrix.transpose m) (Matrix.transpose b)
   in  approxMatrix (selectReal 1e-1 1e-5)
         (Vector.toComplex a)
         (solveLeft (vr <#> diagonal d) vr)



checkForAll ::
   (Show a, QC.Testable test) =>
   Gen.T tag dim a -> (a -> test) -> Tagged tag QC.Property
checkForAll gen = Util.checkForAll (Gen.run gen 3 5)


testsVar ::
   (Show a, Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   ("multiplySquare",
      checkForAll Gen.square multiplySquare) :
   ("multiplyPower",
      Gen.withExtra checkForAll (QC.choose (0,10)) Gen.square multiplyPower) :
   ("multiplyInverse",
      checkForAll Gen.invertible multiplyInverse) :
   ("determinantSingleton",
      checkForAll Gen.scalar determinantSingleton) :
   ("determinantTranspose",
      checkForAll Gen.square determinantTranspose) :
   ("multiplyDeterminant",
      checkForAll ((,) <$> Gen.square <|*|> Gen.square) multiplyDeterminant) :
   ("multiplySolve",
      checkForAll ((,) <$> Gen.invertible <|\|> Gen.matrix) multiplySolve) :

   ("schur",
      checkForAll Gen.square schur) :
   ("eigensystem",
      checkForAll genDiagonalizable eigensystem) :
   ("eigensystemLeft",
      checkForAll genDiagonalizable eigensystemLeft) :
   ("eigensystemRight",
      checkForAll genDiagonalizable eigensystemRight) :
   []
