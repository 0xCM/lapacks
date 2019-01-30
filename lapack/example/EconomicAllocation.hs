{-# LANGUAGE TypeOperators #-}
module Main where

import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix (ZeroInt, (#>), (|||))
import Numeric.LAPACK.Format ((##))

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Shape ((:+:)((:+:)))
import Data.Function.HT (nest)


type ZeroInt2 = ZeroInt:+:ZeroInt
type Vector sh = Vector.Vector sh Double
type Matrix height width = Matrix.General height width Double
type SquareMatrix size = Square.Square size Double


balances :: Vector ZeroInt2
balances =
   Vector.fromList (Matrix.zeroInt 2 :+: Matrix.zeroInt 2)
      [100000, 90000, -50000, -120000]

expenses :: Matrix ZeroInt ZeroInt2
expenses =
   Matrix.fromList (Matrix.zeroInt 2) (Matrix.zeroInt 2 :+: Matrix.zeroInt 2) $
   [16000,  4000,  8000, 12000,
    10000, 30000, 40000, 20000]

normalize ::
   (Eq height, Shape.C height, Shape.C width) =>
   Matrix height width -> Matrix height width
normalize x = Matrix.scaleRows (Array.map recip (Matrix.rowSums x)) x


subtractIdentity :: (Eq sh, Shape.C sh) => SquareMatrix sh -> SquareMatrix sh
subtractIdentity x = Matrix.sub x $ Square.identityFrom x

completeIdSquare :: Matrix ZeroInt2 ZeroInt -> SquareMatrix ZeroInt2
completeIdSquare x =
   Square.fromGeneral $
      (Matrix.takeLeftColumns $ Matrix.fromFull $ Square.identityFromHeight x)
      |||
      x

iterationMatrix :: SquareMatrix ZeroInt2
iterationMatrix =
   completeIdSquare $ Matrix.transpose $ normalize expenses

iterated :: Vector ZeroInt2
iterated = nest 30 (iterationMatrix #>) balances



compensated :: Vector ZeroInt
compensated =
   let a = Matrix.transpose $ normalize expenses
       p = Matrix.takeTopRows a
       k = Square.fromGeneral $ Matrix.takeBottomRows a
       x = Vector.takeLeft balances
       y = Vector.takeRight balances
   in Vector.sub x $ p #> Matrix.solveVector (subtractIdentity k) y


main :: IO ()
main = do
   Array.mapShape (Shape.ZeroBased . Shape.size) iterated ## "%10.2f"
   compensated ## "%10.2f"
