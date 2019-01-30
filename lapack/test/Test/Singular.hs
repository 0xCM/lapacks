{-# LANGUAGE TypeFamilies #-}
module Test.Singular (testsVar) where

import qualified Test.Generator as Gen
import qualified Test.Utility as Util
import Test.Generator ((<|\|>))
import Test.Utility
         (approxReal, approxArrayTol, approxMatrix, isIdentity, Tagged)

import qualified Numeric.LAPACK.Singular as Singular
import qualified Numeric.LAPACK.Orthogonal as Ortho
import qualified Numeric.LAPACK.Matrix.Hermitian as Herm
import qualified Numeric.LAPACK.Matrix as Matrix
import Numeric.LAPACK.Matrix (General, ZeroInt, (<#>))
import Numeric.LAPACK.Scalar (RealOf, selectReal)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape

import Control.Applicative ((<$>))

import qualified Test.QuickCheck as QC


pseudoInverseOrtho ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General ZeroInt ZeroInt a -> Bool
pseudoInverseOrtho a =
   let (no,invo) = Ortho.pseudoInverseRCond 1e-5 a
       (ns,invs) = Singular.pseudoInverseRCond 1e-5 a
       tol = selectReal 1e-2 1e-5
   in no==ns && approxMatrix tol invo invs

pseudoInverseProjection ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General ZeroInt ZeroInt a -> Bool
pseudoInverseProjection a =
   let ainv = snd $ Singular.pseudoInverseRCond 1e-5 a
       tol = selectReal 1e-1 1e-5
   in approxArrayTol tol a (a <#> ainv <#> a) &&
      approxArrayTol tol ainv (ainv <#> a <#> ainv)

pseudoInverseHermitian ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General ZeroInt ZeroInt a -> Bool
pseudoInverseHermitian a =
   let ainv = snd $ Singular.pseudoInverseRCond 1e-5 a
       tol = selectReal 1e-2 1e-5
       aainv = a <#> ainv
       ainva = ainv <#> a
   in approxMatrix tol aainv (Matrix.adjoint aainv) &&
      approxMatrix tol ainva (Matrix.adjoint ainva)


determinantAbsolute ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   General ZeroInt ZeroInt a -> Bool
determinantAbsolute a =
   let detOrtho = Ortho.determinantAbsolute a
       detSing = Singular.determinantAbsolute a
   in approxReal
         (selectReal 1e-3 1e-5 * max 1 (max detOrtho detSing))
         detOrtho detSing


leastSquaresMinimumNorm ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Matrix.General ZeroInt ZeroInt a, Matrix.General ZeroInt ZeroInt a) -> Bool
leastSquaresMinimumNorm (a,b) =
   let (no,xo) = Ortho.leastSquaresMinimumNormRCond 1e-5 a b
       (ns,xs) = Singular.leastSquaresMinimumNormRCond 1e-5 a b
   in no==ns &&
      approxMatrix (selectReal 10 1e-3) xo xs


decompose ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.General ZeroInt ZeroInt a -> Bool
decompose a =
   let (u,s,vt) = Singular.decompose a
       mn = Shape.size $ Array.shape s
   in approxArrayTol 1e-3 a
        (Matrix.takeColumns mn (Matrix.generalizeWide u) <#>
         Matrix.scaleRowsReal s (Matrix.takeRows mn (Matrix.generalizeTall vt)))
      &&
      isIdentity 1e-3 (Matrix.adjoint u <#> u)
      &&
      isIdentity 1e-3 (Matrix.adjoint vt <#> vt)

decomposeTall ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.Tall ZeroInt ZeroInt a -> Bool
decomposeTall a =
   let (u,s,vt) = Singular.decomposeTall a
   in approxArrayTol 1e-3 a (u <#> Matrix.scaleRowsReal s vt)
      &&
      isIdentity 1e-3 (Herm.toSquare $ Herm.covariance $ Matrix.fromFull u)
      &&
      isIdentity 1e-3 (Matrix.adjoint vt <#> vt)

decomposeWide ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   Matrix.Wide ZeroInt ZeroInt a -> Bool
decomposeWide a =
   let (u,s,vt) = Singular.decomposeWide a
   in approxArrayTol 1e-3 a (u <#> Matrix.scaleRowsReal s vt)
      &&
      isIdentity 1e-3 (Matrix.adjoint u <#> u)
      &&
      isIdentity 1e-3
         (Herm.toSquare $ Herm.covariance $
          Matrix.fromFull $ Matrix.transpose vt)



checkForAll ::
   (Show a, QC.Testable test) =>
   Gen.T tag dim a -> (a -> test) -> Tagged tag QC.Property
checkForAll gen = Util.checkForAll (Gen.run gen 3 5)

testsVar ::
   (Show a, Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   ("pseudoInverseOrtho",
      checkForAll Gen.matrix pseudoInverseOrtho) :
   ("pseudoInverseProjection",
      checkForAll Gen.matrix pseudoInverseProjection) :
   ("pseudoInverseHermitian",
      checkForAll Gen.matrix pseudoInverseHermitian) :
   ("determinantAbsolute",
      checkForAll Gen.matrix determinantAbsolute) :
   ("leastSquaresMinimumNorm",
      checkForAll ((,) <$> Gen.matrix <|\|> Gen.matrix) leastSquaresMinimumNorm) :
   ("decompose",
      checkForAll Gen.matrix decompose) :
   ("decomposeTall",
      checkForAll Gen.tall decomposeTall) :
   ("decomposeWide",
      checkForAll Gen.wide decomposeWide) :
   []
