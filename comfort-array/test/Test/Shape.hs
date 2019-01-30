module Test.Shape where

import qualified Data.Array.Comfort.Shape.Test as ShapeTest
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Shape ((:+:)((:+:)))

import qualified Test.QuickCheck as QC
import Test.Utility (prefix)

import Control.Applicative (liftA2, liftA3)


genZeroBased :: Int -> QC.Gen (Shape.ZeroBased Int)
genZeroBased n = fmap Shape.ZeroBased $ QC.choose (0,n)

tests :: [(String, QC.Property)]
tests =
   prefix "ZeroBased"
      (ShapeTest.tests $ genZeroBased 10) ++
   prefix "OneBased"
      (ShapeTest.tests $ fmap Shape.OneBased $ QC.choose (0,10::Int)) ++
   prefix "Range"
      (ShapeTest.tests $ do
         from <- QC.choose (0,10::Int)
         to <- QC.choose (from-1, 10)
         return $ Shape.Range from to) ++
   prefix "Shifted"
      (ShapeTest.tests $
       liftA2 Shape.Shifted
         (QC.choose (-10,10::Int)) (QC.choose (0,10::Int))) ++
   prefix "Enumeration Bool"
      (ShapeTest.tests $
       return (Shape.Enumeration :: Shape.Enumeration Bool)) ++
   prefix "Enumeration Ordering"
      (ShapeTest.tests $
       return (Shape.Enumeration :: Shape.Enumeration Ordering)) ++
   prefix "Deferred Shifted"
      (ShapeTest.tests $ fmap Shape.Deferred $
       liftA2 Shape.Shifted
         (QC.choose (-10,10::Int)) (QC.choose (0,10::Int))) ++
   prefix "Pair"
      (ShapeTest.tests $ liftA2 (,) (genZeroBased 10) (genZeroBased 10)) ++
   prefix "Triple"
      (ShapeTest.tests $
       liftA3 (,,) (genZeroBased 10) (genZeroBased 10) (genZeroBased 10)) ++
   prefix "Append"
      (ShapeTest.tests $ liftA2 (:+:) (genZeroBased 10) (genZeroBased 10)) ++
   prefix "Triangular Lower"
      (ShapeTest.tests $
       fmap (Shape.Triangular Shape.Lower) (genZeroBased 10)) ++
   prefix "Triangular Upper"
      (ShapeTest.tests $
       fmap (Shape.Triangular Shape.Upper) (genZeroBased 10)) ++
   []
