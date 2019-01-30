{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
module Test.Matrix (testsVar) where

import qualified Test.Generator as Gen
import qualified Test.Utility as Util
import Test.Generator
         ((<|*|>), (<|*.>), (<.*.>), (<***>), (<|=|>), (<|||>), (<===>))
import Test.Utility
         (approx, approxArray, approxMatrix,
          genOrder, Tagged(Tagged), TaggedGen)

import qualified Numeric.LAPACK.Matrix.Triangular as Triangular
import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix.Square as Square
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import Numeric.LAPACK.Matrix
         (General, ZeroInt, zeroInt, (#>), (<#>), (|||), (===))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (RealOf, conjugate)

import qualified Numeric.Netlib.Class as Class

import qualified Data.Array.Comfort.Storable as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)
import Data.Array.Comfort.Shape ((:+:))

import Control.Applicative (liftA2, (<$>))

import Data.Tuple.HT (mapPair, swap)
import Data.Eq.HT (equating)

import qualified Test.QuickCheck as QC


genArray ::
   (Shape.C shape, Class.Floating a) => shape -> QC.Gen (Array shape a)
genArray = Util.genArray 10

equalArray ::
   (Shape.C shape, Eq shape, Class.Floating a) =>
   Array shape a -> Array shape a -> Bool
equalArray x y =
   if Array.shape x == Array.shape y
     then equalArrayBody x y
     else error "equalArray: shapes mismatch"

equalArrayBody ::
   (Shape.C shape, Class.Floating a) =>
   Array shape a -> Array shape a -> Bool
equalArrayBody =
   getEqualArray $
   Class.switchFloating
      (EqualArray $ equating Array.toList)
      (EqualArray $ equating Array.toList)
      (EqualArray $ equating Array.toList)
      (EqualArray $ equating Array.toList)

newtype EqualArray f a = EqualArray {getEqualArray :: f a -> f a -> Bool}


dotProduct ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Vector ZeroInt a, Vector ZeroInt a) -> Bool
dotProduct (x,y) =
   approx 1e-5
      (Vector.dot x y)
      (Matrix.toScalar $
       Matrix.singleRow MatrixShape.RowMajor x <#>
       Matrix.singleColumn MatrixShape.ColumnMajor y)

innerDot ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Vector ZeroInt a, Vector ZeroInt a) -> Bool
innerDot (x,y) =
   approx 1e-5 (Vector.inner x y) (Vector.dot (Vector.conjugate x) y)

tensorProductTranspose ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order -> (Vector ZeroInt a, Vector ZeroInt a) -> Bool
tensorProductTranspose order (x,y) =
   approxArray
      (Matrix.transpose (Matrix.tensorProduct order x y))
      (Matrix.tensorProduct (MatrixShape.flipOrder order) y x)

outerTranspose ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order -> (Vector ZeroInt a, Vector ZeroInt a) -> Bool
outerTranspose order (x,y) =
   approxArray
      (Matrix.transpose (Matrix.outer order x y))
      (Matrix.outer (MatrixShape.flipOrder order)
         (Vector.conjugate y) (Vector.conjugate x))

tensorProduct ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order -> (Vector ZeroInt a, Vector ZeroInt a) -> Bool
tensorProduct order (x,y) =
   approxArray
      (Matrix.tensorProduct order x y)
      (Matrix.singleColumn order x <#> Matrix.singleRow order y)

tensorProductMul ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular.Diagonal ZeroInt a,
    Matrix.General ZeroInt ZeroInt a,
    Triangular.Diagonal ZeroInt a) ->
   Bool
tensorProductMul (x,m,y) =
   let xmy = x <#> m <#> y
   in approxArray xmy
         (Vector.mul m
            (Matrix.tensorProduct (MatrixShape.fullOrder $ Array.shape xmy)
               (Triangular.takeDiagonal x) (Triangular.takeDiagonal y)))

outerTensorProduct ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order -> (Vector ZeroInt a, Vector ZeroInt a) -> Bool
outerTensorProduct order (x,y) =
   approxArray
      (Matrix.outer order x y)
      (Matrix.tensorProduct order x $ Vector.conjugate y)

genScaledVectorPairs ::
   (Class.Floating a) =>
   Gen.Matrix a Int Int
      ((ZeroInt, ZeroInt), [(a, (Vector ZeroInt a, Vector ZeroInt a))])
genScaledVectorPairs =
   flip Gen.mapGen Gen.matrixDims $ \maxElem size@(height,width) ->
      fmap ((,) size) $
      QC.listOf $
         liftA2 (,) (Util.genElement maxElem) $
         liftA2 (,) (Util.genArray maxElem height) (Util.genArray maxElem width)

sumRank1 ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order ->
   ((ZeroInt,ZeroInt), [(a, (Vector ZeroInt a, Vector ZeroInt a))]) -> Bool
sumRank1 order (size,xys) =
   approxArray
      (case order of
         MatrixShape.ColumnMajor -> Matrix.sumRank1 size xys
         MatrixShape.RowMajor ->
            Matrix.adjoint $
            Matrix.sumRank1 (swap size) $ map (mapPair (conjugate, swap)) xys)
      (foldl Vector.add
         (Vector.constant (uncurry (MatrixShape.general order) size) 0)
         (map (\(a,(x,y)) -> Matrix.outer order (Vector.scale a x) y) xys))


outerTrace ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order -> (Vector ZeroInt a, Vector ZeroInt a) -> Bool
outerTrace order (x,y) =
   approx 1e-5
      (Vector.inner y x)
      (Square.trace $ Square.fromGeneral $ Matrix.outer order x y)

outerInner ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order ->
   (Vector ZeroInt a, Vector ZeroInt a, Vector ZeroInt a) -> Bool
outerInner order (x,y,z) =
   approxArray (Matrix.outer order x y #> z) (Vector.scale (Vector.inner y z) x)


tensorTrace ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order -> (Vector ZeroInt a, Vector ZeroInt a) -> Bool
tensorTrace order (x,y) =
   approx 1e-5 (Vector.dot y x)
      (Square.trace $ Square.fromGeneral $ Matrix.tensorProduct order x y)

tensorDot ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order ->
   (Vector ZeroInt a, Vector ZeroInt a, Vector ZeroInt a) -> Bool
tensorDot order (x,y,z) =
   approxArray
      (Matrix.tensorProduct order x y #> z) (Vector.scale (Vector.dot y z) x)


genZeroColumns ::
   (Class.Floating a) => TaggedGen a (Matrix.Tall ZeroInt ZeroInt a)
genZeroColumns = Tagged $ do
   height <- zeroInt <$> QC.choose (0,5)
   order <- genOrder
   genArray (MatrixShape.tall order height (zeroInt 0))


reverseNoRows :: (Class.Floating a) => Matrix.Wide ZeroInt ZeroInt a -> Bool
reverseNoRows x =
   equalArray x $ Matrix.reverseRows x

reverseNoColumns :: (Class.Floating a) => Matrix.Tall ZeroInt ZeroInt a -> Bool
reverseNoColumns x =
   equalArray x $ Matrix.reverseColumns x



genMatrix2EqHeight ::
   (Class.Floating a) =>
   Gen.Matrix a Int (Int:+:Int)
      (General ZeroInt ZeroInt a, General ZeroInt ZeroInt a)
genMatrix2EqHeight = (,) <$> Gen.matrix <|||> Gen.matrix

genMatrix2EqWidth ::
   (Class.Floating a) =>
   Gen.Matrix a (Int:+:Int) Int
      (General ZeroInt ZeroInt a, General ZeroInt ZeroInt a)
genMatrix2EqWidth = (,) <$> Gen.matrix <===> Gen.matrix

reverseRows :: (Class.Floating a) => General ZeroInt ZeroInt a -> Bool
reverseRows x =
   equalArray x $ Matrix.reverseRows (Matrix.reverseRows x)

reverseColumns :: (Class.Floating a) => General ZeroInt ZeroInt a -> Bool
reverseColumns x =
   equalArray x $ Matrix.reverseColumns (Matrix.reverseColumns x)


mapHeight ::
   (heightA -> heightB) ->
   MatrixShape.General heightA width ->
   MatrixShape.General heightB width
mapHeight f shape =
   MatrixShape.general
      (MatrixShape.fullOrder shape)
      (f $ MatrixShape.fullHeight shape)
      (MatrixShape.fullWidth shape)

mapWidth ::
   (widthA -> widthB) ->
   MatrixShape.General height widthA ->
   MatrixShape.General height widthB
mapWidth f shape =
   MatrixShape.general
      (MatrixShape.fullOrder shape)
      (MatrixShape.fullHeight shape)
      (f $ MatrixShape.fullWidth shape)

zeroIntHeight ::
   (Shape.C height, Shape.C width) =>
   General height width a -> General ZeroInt width a
zeroIntHeight = Array.mapShape (mapHeight (zeroInt . Shape.size))

zeroIntWidth ::
   (Shape.C height, Shape.C width) =>
   General height width a -> General height ZeroInt a
zeroIntWidth = Array.mapShape (mapWidth (zeroInt . Shape.size))

reverseRowsStack ::
   (Class.Floating a) =>
   (General ZeroInt ZeroInt a, General ZeroInt ZeroInt a) -> Bool
reverseRowsStack (x,y) =
   equalArray
      (Matrix.reverseRows $ zeroIntHeight $ x===y)
      (zeroIntHeight $ Matrix.reverseRows y === Matrix.reverseRows x)

reverseColumnsStack ::
   (Class.Floating a) =>
   (General ZeroInt ZeroInt a, General ZeroInt ZeroInt a) -> Bool
reverseColumnsStack (x,y) =
   equalArray
      (Matrix.reverseColumns $ zeroIntWidth $ x|||y)
      (zeroIntWidth $ Matrix.reverseColumns y ||| Matrix.reverseColumns x)


data Cut = Take | Drop deriving (Show, Eq, Ord, Enum, Bounded)
data Slice = Row | Column deriving (Show, Eq, Ord, Enum, Bounded)

cut ::
   (Class.Floating a) =>
   Cut -> Slice -> Int ->
   General ZeroInt ZeroInt a -> General ZeroInt ZeroInt a
cut Take Row = Matrix.takeRows
cut Take Column = Matrix.takeColumns
cut Drop Row = Matrix.dropRows
cut Drop Column = Matrix.dropColumns

cutCommutative ::
   (Class.Floating a) =>
   ((Cut,Slice),(Int,Int)) -> General ZeroInt ZeroInt a -> Bool
cutCommutative (kind,(k,j)) x =
   let cutK = uncurry cut kind k
       cutJ = uncurry cut kind j
   in equalArray (cutK $ cutJ x) (cutJ $ cutK x)

cutRowColumnCommutative ::
   (Class.Floating a) =>
   ((Cut,Int),(Cut,Int)) -> General ZeroInt ZeroInt a -> Bool
cutRowColumnCommutative ((cutR,k),(cutC,j)) x =
   let cutRows = cut cutR Row k
       cutColumns = cut cutC Column j
   in equalArray (cutRows $ cutColumns x) (cutColumns $ cutRows x)


takeEqually ::
   (Class.Floating a) => Int -> General ZeroInt ZeroInt a -> Bool
takeEqually k x =
   equalArray
      (Matrix.takeEqually k x)
      (Matrix.takeRows k (Matrix.takeColumns k x))

dropEqually ::
   (Class.Floating a) => Int -> General ZeroInt ZeroInt a -> Bool
dropEqually k x =
   equalArray
      (Matrix.dropEqually k x)
      (Matrix.dropRows k (Matrix.dropColumns k x))


stackSplitRows ::
   (Class.Floating a) => Int -> General ZeroInt ZeroInt a -> Bool
stackSplitRows k x =
   equalArray x
      (zeroIntHeight $ Matrix.takeRows k x === Matrix.dropRows k x)

stackSplitColumns ::
   (Class.Floating a) => Int -> General ZeroInt ZeroInt a -> Bool
stackSplitColumns k x =
   equalArray x
      (zeroIntWidth $ Matrix.takeColumns k x ||| Matrix.dropColumns k x)


takeStackRows, dropStackRows ::
   (Class.Floating a) =>
   (General ZeroInt ZeroInt a, General ZeroInt ZeroInt a) -> Bool
takeStackRows (x,y) =
   equalArray
      (Matrix.toRowMajor x)
      (Matrix.toRowMajor $ Matrix.takeRows (Shape.size $ Matrix.height x) $
       zeroIntHeight $ x===y)
dropStackRows (x,y) =
   equalArray
      (Matrix.toRowMajor y)
      (Matrix.toRowMajor $ Matrix.dropRows (Shape.size $ Matrix.height x) $
       zeroIntHeight $ x===y)

takeStackColumns, dropStackColumns ::
   (Class.Floating a) =>
   (General ZeroInt ZeroInt a, General ZeroInt ZeroInt a) -> Bool
takeStackColumns (x,y) =
   equalArray
      (Matrix.toRowMajor x)
      (Matrix.toRowMajor $ Matrix.takeColumns (Shape.size $ Matrix.width x) $
       zeroIntWidth $ x|||y)
dropStackColumns (x,y) =
   equalArray
      (Matrix.toRowMajor y)
      (Matrix.toRowMajor $ Matrix.dropColumns (Shape.size $ Matrix.width x) $
       zeroIntWidth $ x|||y)

stackRowsAssociative, stackColumnsAssociative ::
   (Class.Floating a) =>
   (General ZeroInt ZeroInt a,
    General ZeroInt ZeroInt a,
    General ZeroInt ZeroInt a) -> Bool
stackRowsAssociative (x,y,z) =
   equalArray
      (zeroIntHeight ((x===y)===z))
      (zeroIntHeight (x===(y===z)))
stackColumnsAssociative (x,y,z) =
   equalArray
      (zeroIntWidth ((x|||y)|||z))
      (zeroIntWidth (x|||(y|||z)))

stackRowsColumnsCommutative ::
   (Class.Floating a) =>
   ((General ZeroInt ZeroInt a, General ZeroInt ZeroInt a),
    (General ZeroInt ZeroInt a, General ZeroInt ZeroInt a)) -> Bool
stackRowsColumnsCommutative ((x,y),(z,w)) =
   equalArray
      (Matrix.toRowMajor $ (x|||y)===(z|||w))
      (Matrix.toRowMajor $ (x===z)|||(y===w))


forceOrder ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   MatrixShape.Order ->
   (General ZeroInt ZeroInt a, Vector ZeroInt a) ->
   Bool
forceOrder order (a,x) =
   let ao = Matrix.forceOrder order a
   in MatrixShape.fullOrder (Array.shape ao) == order
      &&
      approxArray (a #> x) (ao #> x)

addDistributive ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   ((General ZeroInt ZeroInt a, General ZeroInt ZeroInt a), Vector ZeroInt a) ->
   Bool
addDistributive ((a,b),x) =
   approxArray (Matrix.add a b #> x) (Vector.add (a#>x) (b#>x))

subDistributive ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   ((General ZeroInt ZeroInt a, General ZeroInt ZeroInt a), Vector ZeroInt a) ->
   Bool
subDistributive ((a,b),x) =
   approxArray (Matrix.sub a b #> x) (Vector.sub (a#>x) (b#>x))


multiplyDiagonalMatrix ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (Triangular.Diagonal ZeroInt a, General ZeroInt ZeroInt a) -> Bool
multiplyDiagonalMatrix (x,y) =
   approxArray (x <#> y) (Triangular.toSquare x <#> y)

multiplyMatrixDiagonal ::
   (Class.Floating a, RealOf a ~ ar, Class.Real ar) =>
   (General ZeroInt ZeroInt a, Triangular.Diagonal ZeroInt a) -> Bool
multiplyMatrixDiagonal (x,y) =
   approxMatrix 1e-5 (x <#> y) (x <#> Triangular.toSquare y)



checkForAll ::
   (Show a, QC.Testable test) =>
   Gen.T tag dim a -> (a -> test) -> Tagged tag QC.Property
checkForAll gen = Util.checkForAll (Gen.run gen 10 5)

checkForAllExtra ::
   (Show a, Show b, QC.Testable test) =>
   QC.Gen a -> Gen.T tag dim b ->
   (a -> b -> test) -> Tagged tag QC.Property
checkForAllExtra = Gen.withExtra checkForAll


testsVar ::
   (Show a, Class.Floating a, Eq a, RealOf a ~ ar, Class.Real ar) =>
   [(String, Tagged a QC.Property)]
testsVar =
   ("dotProduct",
      checkForAll ((,) <$> Gen.vector <.*.> Gen.vector) dotProduct) :
   ("innerDot",
      checkForAll ((,) <$> Gen.vector <.*.> Gen.vector) innerDot) :
   ("tensorProductTranspose",
      checkForAllExtra genOrder
         ((,) <$> Gen.vector <***> Gen.vector) tensorProductTranspose) :
   ("outerTranspose",
      checkForAllExtra genOrder
         ((,) <$> Gen.vector <***> Gen.vector) outerTranspose) :
   ("tensorProduct",
      checkForAllExtra genOrder
         ((,) <$> Gen.vector <***> Gen.vector) tensorProduct) :
   ("tensorProductMul",
      checkForAll ((,,) <$> Gen.diagonal <|*|> Gen.matrix <|*|> Gen.diagonal)
         tensorProductMul) :
   ("outerTensorProduct",
      checkForAllExtra genOrder
         ((,) <$> Gen.vector <***> Gen.vector) outerTensorProduct) :
   ("sumRank1",
      checkForAllExtra genOrder genScaledVectorPairs sumRank1) :

   ("outerTrace",
      checkForAllExtra genOrder
         ((,) <$> Gen.vector <.*.> Gen.vector) outerTrace) :
   ("outerInner",
      checkForAllExtra genOrder
         ((,,) <$> Gen.vector <***> Gen.vector <|*.> Gen.vector) outerInner) :
   ("tensorTrace",
      checkForAllExtra genOrder
         ((,) <$> Gen.vector <.*.> Gen.vector) tensorTrace) :
   ("tensorDot",
      checkForAllExtra genOrder
         ((,,) <$> Gen.vector <***> Gen.vector <|*.> Gen.vector) tensorDot) :

   ("reverseNoRows",
      Util.checkForAllPlain
         (fmap Matrix.transpose <$> genZeroColumns) reverseNoRows) :
   ("reverseNoColumns",
      Util.checkForAllPlain genZeroColumns reverseNoColumns) :
   ("reverseRows",
      checkForAll Gen.matrix reverseRows) :
   ("reverseColumns",
      checkForAll Gen.matrix reverseColumns) :
   ("reverseRowsStack",
      checkForAll genMatrix2EqWidth reverseRowsStack) :
   ("reverseColumnsStack",
      checkForAll genMatrix2EqHeight reverseColumnsStack) :
   ("cutCommutative",
      checkForAllExtra
         (liftA2 (,)
            (liftA2 (,) QC.arbitraryBoundedEnum QC.arbitraryBoundedEnum)
            (liftA2 (,) (QC.choose (0,5)) (QC.choose (0,5))))
         Gen.matrix cutCommutative) :
   ("cutRowColumnCommutative",
      checkForAllExtra
         (liftA2 (,)
            (liftA2 (,) QC.arbitraryBoundedEnum (QC.choose (0,5)))
            (liftA2 (,) QC.arbitraryBoundedEnum (QC.choose (0,5))))
         Gen.matrix cutRowColumnCommutative) :
   ("takeEqually",
      checkForAllExtra (QC.choose (0,5)) Gen.matrix takeEqually) :
   ("dropEqually",
      checkForAllExtra (QC.choose (0,5)) Gen.matrix dropEqually) :
   ("stackSplitRows",
      checkForAllExtra (QC.choose (0,5)) Gen.matrix stackSplitRows) :
   ("stackSplitColumns",
      checkForAllExtra (QC.choose (0,5)) Gen.matrix stackSplitColumns) :
   ("takeStackRows",
      checkForAll genMatrix2EqWidth takeStackRows) :
   ("dropStackRows",
      checkForAll genMatrix2EqWidth dropStackRows) :
   ("takeStackColumns",
      checkForAll genMatrix2EqHeight takeStackColumns) :
   ("dropStackColumns",
      checkForAll genMatrix2EqHeight dropStackColumns) :
   ("stackRowsAssociative",
      checkForAll
         ((,,) <$> Gen.matrix <===> Gen.matrix <===> Gen.matrix)
         stackRowsAssociative) :
   ("stackColumnsAssociative",
      checkForAll
         ((,,) <$> Gen.matrix <|||> Gen.matrix <|||> Gen.matrix)
         stackColumnsAssociative) :
   ("stackRowsColumnsCommutative",
      checkForAll
         ((,) <$> genMatrix2EqHeight <===> genMatrix2EqHeight)
         stackRowsColumnsCommutative) :

   ("forceOrder",
      checkForAllExtra genOrder
         ((,) <$> Gen.matrix <|*.> Gen.vector) forceOrder) :
   ("addDistributive",
      checkForAll
         ((,) <$>
            ((,) <$> Gen.matrix <|=|> Gen.matrix)
            <|*.>
            Gen.vector)
         addDistributive) :
   ("subDistributive",
      checkForAll
         ((,) <$>
            ((,) <$> Gen.matrix <|=|> Gen.matrix)
            <|*.>
            Gen.vector)
         subDistributive) :

   ("multiplyDiagonalMatrix",
      checkForAll
         ((,) <$> Gen.diagonal <|*|> Gen.matrix) multiplyDiagonalMatrix) :
   ("multiplyMatrixDiagonal",
      checkForAll
         ((,) <$> Gen.matrix <|*|> Gen.diagonal) multiplyMatrixDiagonal) :
   []
