{-# LANGUAGE Rank2Types #-}
module Test.Format where

import qualified Numeric.LAPACK.Matrix.Shape as MatrixShape
import qualified Numeric.LAPACK.Matrix.BandedHermitian as BandedHermitian
import qualified Numeric.LAPACK.Matrix.Banded as Banded
import qualified Numeric.LAPACK.Matrix.Hermitian as Hermitian
import qualified Numeric.LAPACK.Matrix.Triangular as Triangular
import qualified Numeric.LAPACK.Matrix as Matrix
import qualified Numeric.LAPACK.Vector as Vector
import qualified Numeric.LAPACK.Permutation as Perm
import Numeric.LAPACK.Matrix.Shape (Order(RowMajor, ColumnMajor), UnaryProxy)
import Numeric.LAPACK.Matrix (ZeroInt, zeroInt)
import Numeric.LAPACK.Format (Format, (##))

import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary.Literal as TypeNum
import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num.Unary (unary)

import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)

import Data.Complex as Cplx (Complex((:+)))


vector :: (Class.Floating a) => Vector.Vector ZeroInt a
vector = Vector.random Vector.UniformBoxPM1 (zeroInt 4) 419

general :: (Class.Floating a) => Order -> Matrix.General ZeroInt ZeroInt a
general order =
   Vector.random Vector.UniformBoxPM1
      (MatrixShape.general order (zeroInt 3) (zeroInt 4)) 420

split ::
   (Eq lower, Shape.C height, Shape.C width, Class.Floating a) =>
   lower -> height -> width -> Order ->
   Array (MatrixShape.SplitGeneral lower height width) a
split lowerPart height width order =
   Vector.random Vector.UniformBoxPM1
      (MatrixShape.splitGeneral lowerPart order height width) 420

hermitian :: (Class.Floating a) => Order -> Hermitian.Hermitian ZeroInt a
hermitian order =
   Vector.random Vector.UniformBoxPM1
      (MatrixShape.hermitian order (zeroInt 4)) 421

diagonal :: (Class.Floating a) => Order -> Triangular.Diagonal ZeroInt a
diagonal order =
   Vector.random Vector.UniformBoxPM1
      (MatrixShape.diagonal order (zeroInt 4)) 422

lowerTriangular ::
   (Class.Floating a) => Order -> Triangular.Lower ZeroInt a
lowerTriangular order =
   Vector.random Vector.UniformBoxPM1
      (MatrixShape.lowerTriangular order (zeroInt 4)) 423

upperTriangular ::
   (Class.Floating a) => Order -> Triangular.Upper ZeroInt a
upperTriangular order =
   Vector.random Vector.UniformBoxPM1
      (MatrixShape.upperTriangular order (zeroInt 4)) 424

symmetric :: (Class.Floating a) => Order -> Triangular.Symmetric ZeroInt a
symmetric order =
   Vector.random Vector.UniformBoxPM1
      (MatrixShape.symmetric order (zeroInt 4)) 425


bandedHermitian ::
   (Unary.Natural offDiag, Class.Floating a) =>
   UnaryProxy offDiag -> Order ->
   BandedHermitian.BandedHermitian offDiag ZeroInt a
bandedHermitian numOff order =
   Vector.random Vector.UniformBoxPM1
      (MatrixShape.bandedHermitian numOff order (zeroInt 4)) 426

banded ::
   (Unary.Natural sub, Unary.Natural super,
    Shape.C height, Shape.C width, Class.Floating a) =>
   (UnaryProxy sub, UnaryProxy super) -> height -> width -> Order ->
   Banded.General sub super height width a
banded offDiag height width order =
   Vector.random Vector.UniformBoxPM1
      (MatrixShape.bandedGeneral offDiag order height width) 427


permutation :: Perm.Permutation ZeroInt
permutation =
   Perm.fromPivots Perm.NonInverted (zeroInt 5) $
   Vector.fromList (zeroInt 5) [3,2,4,5,5]


fmt :: String
fmt = "%.4g"

printFormatted :: Format a => a -> IO ()
printFormatted x = putStrLn "" >> (x ## fmt)

printVectorFloat :: (Format (f Float)) => f Float -> IO ()
printVectorFloat = printFormatted

printVectorComplex ::
   (Format (f (Complex Float))) => f (Complex Float) -> IO ()
printVectorComplex = printFormatted

printVectorWithOrder ::
   Format (f Float) =>
   Format (f (Complex Float)) =>
   (forall a. (Class.Floating a) => Order -> f a) -> IO ()
printVectorWithOrder f = do
   printFormatted $ floatVector $ f RowMajor
   printFormatted $ floatVector $ f ColumnMajor
   printFormatted $ complexVector $ f RowMajor
   printFormatted $ complexVector $ f ColumnMajor

floatVector :: f Float -> f Float
floatVector = id

complexVector :: f (Complex Float) -> f (Complex Float)
complexVector = id

main :: IO ()
main = do
   printFormatted (pi :: Float)
   printFormatted permutation
   printVectorFloat $ sin (1:+1)
   printVectorFloat vector
   printVectorComplex vector
   printVectorWithOrder general
   printVectorWithOrder $ split MatrixShape.Reflector (zeroInt 4) (zeroInt 3)
   printVectorWithOrder $ split MatrixShape.Reflector (zeroInt 3) (zeroInt 4)
   printVectorWithOrder $ split MatrixShape.Triangle (zeroInt 4) (zeroInt 3)
   printVectorWithOrder hermitian
   printVectorWithOrder diagonal
   printVectorWithOrder lowerTriangular
   printVectorWithOrder upperTriangular
   printVectorWithOrder symmetric
   printVectorWithOrder $ bandedHermitian $ unary TypeNum.u0
   printVectorWithOrder $ bandedHermitian $ unary TypeNum.u1
   printVectorWithOrder $ bandedHermitian $ unary TypeNum.u2
   printVectorWithOrder $
      banded (unary TypeNum.u0, unary TypeNum.u0) (zeroInt 4) (zeroInt 3)
   printVectorWithOrder $
      banded (unary TypeNum.u0, unary TypeNum.u2) (zeroInt 4) (zeroInt 3)
   printVectorWithOrder $
      banded (unary TypeNum.u2, unary TypeNum.u0) (zeroInt 4) (zeroInt 3)
   printVectorWithOrder $
      banded (unary TypeNum.u1, unary TypeNum.u2) (zeroInt 4) (zeroInt 3)
   printVectorWithOrder $
      banded (unary TypeNum.u1, unary TypeNum.u2) (zeroInt 3) (zeroInt 4)
