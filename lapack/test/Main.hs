{-# LANGUAGE TypeFamilies #-}
module Main where

import qualified Test.Vector as Vector
import qualified Test.Matrix as Matrix
import qualified Test.Square as Square
import qualified Test.Triangular as Triangular
import qualified Test.Hermitian as Hermitian
import qualified Test.Banded as Banded
import qualified Test.BandedHermitian as BandedHermitian
import qualified Test.Orthogonal as Orthogonal
import qualified Test.Singular as Singular
import qualified Test.Shape as Shape
import qualified Test.Permutation as Permutation
import Test.Format ()
import Test.Utility (Tagged(Tagged), prefix)

import qualified Test.QuickCheck as QC

import Numeric.LAPACK.Scalar (RealOf)

import qualified Numeric.Netlib.Class as Class

import Type.Base.Proxy (Proxy(Proxy))

import qualified Data.List as List
import Data.Complex (Complex)
import Data.Tuple.HT (mapSnd)


testsVar ::
   (Show a, Show ar,
    Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar, RealOf ar ~ ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   prefix "Vector" Vector.testsVar ++
   prefix "Matrix" Matrix.testsVar ++
   prefix "Square" Square.testsVar ++
   prefix "Triangular" Triangular.testsVar ++
   prefix "Hermitian" Hermitian.testsVar ++
   prefix "Banded" Banded.testsVar ++
   prefix "BandedHermitian" BandedHermitian.testsVar ++
   prefix "Orthogonal" Orthogonal.testsVar ++
   prefix "Singular" Singular.testsVar ++
   []

tagTests ::
   String -> Proxy tag ->
   [(String, Tagged tag QC.Property)] -> [(String, QC.Property)]
tagTests typeName Proxy =
   map (\(name, Tagged prop) -> (name++"."++typeName, prop))

tests :: [(String, QC.Property)]
tests =
   concat $ List.transpose $
   (tagTests "Float" (Proxy :: Proxy Float) testsVar) :
   (tagTests "Double" (Proxy :: Proxy Double) testsVar) :
   (tagTests "ComplexFloat" (Proxy :: Proxy (Complex Float)) testsVar) :
   (tagTests "ComplexDouble" (Proxy :: Proxy (Complex Double)) testsVar) :
   []

simpleTests :: [(String, QC.Property)]
simpleTests =
   prefix "Shape" Shape.tests ++
   prefix "Permutation" Permutation.tests ++
   []

main :: IO ()
main =
   mapM_ (\(name,act) -> putStr (name ++ ": ") >> act) $

   map (mapSnd (QC.quickCheckWith (QC.stdArgs {QC.maxSuccess=200}))) tests
   ++
   map (mapSnd QC.quickCheck) simpleTests
