{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE UndecidableInstances #-}
module Numeric.LAPACK.Matrix.Multiply where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as ExtentPriv
import qualified Numeric.LAPACK.Matrix.Extent as Extent
import qualified Numeric.LAPACK.Matrix.BandedHermitian.Basic as BandedHermitian
import qualified Numeric.LAPACK.Matrix.Banded.Basic as Banded
import qualified Numeric.LAPACK.Matrix.Triangular.Basic as Triangular
import qualified Numeric.LAPACK.Matrix.Hermitian.Basic as Hermitian
import qualified Numeric.LAPACK.Private as Private
import Numeric.LAPACK.Matrix.Shape.Private
         (HeightOf, WidthOf, Empty, Filled, Unit, NonUnit,
          Order(RowMajor,ColumnMajor), flipOrder, transposeFromOrder)
import Numeric.LAPACK.Matrix.Extent.Private (Small)
import Numeric.LAPACK.Matrix.Triangular.Basic (Triangular)
import Numeric.LAPACK.Matrix.Basic (transpose)
import Numeric.LAPACK.Matrix.Private
         (Square, Full, mapExtent,
          Transposition(NonTransposed, Transposed))
import Numeric.LAPACK.Vector (Vector)
import Numeric.LAPACK.Scalar (zero, one)

import qualified Numeric.Netlib.Utility as Call
import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num.Unary ((:+:))

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable.Unchecked (Array(Array))

import Foreign.ForeignPtr (withForeignPtr)

import Control.Monad.Trans.Cont (ContT(ContT), evalContT)
import Control.Monad.IO.Class (liftIO)


multiplyVector ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width, Eq width,
    Class.Floating a) =>
   Full vert horiz height width a -> Vector width a -> Vector height a
multiplyVector a x =
   let width = MatrixShape.fullWidth $ Array.shape a
   in if width == Array.shape x
         then multiplyVectorUnchecked a x
         else error "multiplyVector: width shapes mismatch"

multiplyVectorUnchecked ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   Full vert horiz height width a -> Vector width a -> Vector height a
multiplyVectorUnchecked
   (Array shape@(MatrixShape.Full order extent) a) (Array _ x) =
      Array.unsafeCreate (Extent.height extent) $ \yPtr -> do
   let (m,n) = MatrixShape.dimensions shape
   let lda = m
   evalContT $ do
      transPtr <- Call.char $ transposeFromOrder order
      mPtr <- Call.cint m
      nPtr <- Call.cint n
      alphaPtr <- Call.number one
      aPtr <- ContT $ withForeignPtr a
      ldaPtr <- Call.leadingDim lda
      xPtr <- ContT $ withForeignPtr x
      incxPtr <- Call.cint 1
      betaPtr <- Call.number zero
      incyPtr <- Call.cint 1
      liftIO $
         Private.gemv
            transPtr mPtr nPtr alphaPtr aPtr ldaPtr
            xPtr incxPtr betaPtr yPtr incyPtr

{- |
Multiply two matrices with the same dimension constraints.
E.g. you can multiply 'General' and 'General' matrices,
or 'Square' and 'Square' matrices.
It may seem to be overly strict in this respect,
but that design supports type inference the best.
You can lift the restrictions by generalizing operands
with 'Square.toFull', 'Matrix.fromFull',
'Matrix.generalizeTall' or 'Matrix.generalizeWide'.
-}
multiply, multiplyColumnMajor ::
   (Extent.C vert, Extent.C horiz,
    Shape.C height,
    Shape.C fuse, Eq fuse,
    Shape.C width,
    Class.Floating a) =>
   Full vert horiz height fuse a ->
   Full vert horiz fuse width a ->
   Full vert horiz height width a
-- preserve order of the right factor
multiply
   (Array (MatrixShape.Full orderA extentA) a)
   (Array (MatrixShape.Full orderB extentB) b) =
   case Extent.fuse extentA extentB of
      Nothing -> error "multiply: fuse shapes mismatch"
      Just extent ->
         Array.unsafeCreate (MatrixShape.Full orderB extent) $ \cPtr -> do

      let (height,fuse) = Extent.dimensions extentA
      let width = Extent.width extentB
      let m = Shape.size height
      let n = Shape.size width
      let k = Shape.size fuse
      case orderB of
         RowMajor ->
            Private.multiplyMatrix (flipOrder orderB) (flipOrder orderA)
               n k m b a cPtr
         ColumnMajor -> Private.multiplyMatrix orderA orderB m k n a b cPtr

-- always return ColumnMajor
multiplyColumnMajor
   (Array (MatrixShape.Full orderA extentA) a)
   (Array (MatrixShape.Full orderB extentB) b) =
   case Extent.fuse extentA extentB of
      Nothing -> error "multiply: fuse shapes mismatch"
      Just extent ->
         Array.unsafeCreate (MatrixShape.Full ColumnMajor extent) $ \cPtr -> do

      let (height,fuse) = Extent.dimensions extentA
      let width = Extent.width extentB
      let m = Shape.size height
      let n = Shape.size width
      let k = Shape.size fuse
      Private.multiplyMatrix orderA orderB m k n a b cPtr


infixl 7 <#, <#>
infixr 7 #>

class (Shape.C shape) => MultiplyRight shape where
   (#>) ::
      (Class.Floating a) =>
      Array shape a -> Vector (WidthOf shape) a -> Vector (HeightOf shape) a

class (Shape.C shape) => MultiplyLeft shape where
   (<#) ::
      (Class.Floating a) =>
      Vector (HeightOf shape) a -> Array shape a -> Vector (WidthOf shape) a


instance
   (Extent.C vert, Extent.C horiz, Eq width, Shape.C width, Shape.C height) =>
      MultiplyRight (MatrixShape.Full vert horiz height width) where
   (#>) = multiplyVector

instance
   (Extent.C vert, Extent.C horiz, Eq height, Shape.C width, Shape.C height) =>
      MultiplyLeft (MatrixShape.Full vert horiz height width) where
   v <# m = multiplyVector (transpose m) v


instance
   (Eq shape, Shape.C shape) =>
      MultiplyRight (MatrixShape.Hermitian shape) where
   (#>) = Hermitian.multiplyVector NonTransposed

instance
   (Eq shape, Shape.C shape) =>
      MultiplyLeft (MatrixShape.Hermitian shape) where
   (<#) = flip $ Hermitian.multiplyVector Transposed


instance
   (MatrixShape.Content lo, MatrixShape.Content up,
    MatrixShape.TriDiag diag, Eq shape, Shape.C shape) =>
      MultiplyRight (MatrixShape.Triangular lo diag up shape) where
   m #> v = Triangular.multiplyVector m v

instance
   (MatrixShape.Content lo, MatrixShape.Content up,
    MatrixShape.TriDiag diag, Eq shape, Shape.C shape) =>
      MultiplyLeft (MatrixShape.Triangular lo diag up shape) where
   v <# m = Triangular.multiplyVector (Triangular.transpose m) v


instance
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz, Eq width, Shape.C width, Shape.C height) =>
      MultiplyRight (MatrixShape.Banded sub super vert horiz height width) where
   m #> v = Banded.multiplyVector m v

instance
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz, Eq height, Shape.C width, Shape.C height) =>
      MultiplyLeft (MatrixShape.Banded sub super vert horiz height width) where
   v <# m = Banded.multiplyVector (Banded.transpose m) v


instance
   (Unary.Natural offDiag, Shape.C size, Eq size) =>
      MultiplyRight (MatrixShape.BandedHermitian offDiag size) where
   (#>) = BandedHermitian.multiplyVector NonTransposed

instance
   (Unary.Natural offDiag, Shape.C size, Eq size) =>
      MultiplyLeft (MatrixShape.BandedHermitian offDiag size) where
   (<#) = flip $ BandedHermitian.multiplyVector Transposed


{- |
This class allows to multiply two matrices of arbitrary special features
and returns the most special matrix type possible.
At the first glance, this is handy.
At the second glance, this has some problems.
First of all, we may refine the types in future
and then multiplication may return a different, more special type than before.
Second, if you write code with polymorphic matrix types,
then '<#>' may leave you with constraints like
@ExtentPriv.Multiply vert vert ~ vert@.
That constraint is always fulfilled but the compiler cannot infer that.
Because of these problems
you may instead consider using specialised 'multiply' functions
from the various modules for production use.
Btw. 'MultiplyLeft' and 'MultiplyRight' are much less problematic,
because the input and output are always dense vectors.
-}
class (Shape.C shapeA, Shape.C shapeB) => Multiply shapeA shapeB where
   type Multiplied shapeA shapeB
   (<#>) ::
      (Class.Floating a) =>
      Array shapeA a -> Array shapeB a -> Array (Multiplied shapeA shapeB) a

instance
   (Shape.C heightA, Shape.C widthA, Shape.C widthB,
    widthA ~ heightB, Eq heightB,
    Extent.C vertA, Extent.C horizA, Extent.C vertB, Extent.C horizB) =>
      Multiply
         (MatrixShape.Full vertA horizA heightA widthA)
         (MatrixShape.Full vertB horizB heightB widthB) where
   type Multiplied
         (MatrixShape.Full vertA horizA heightA widthA)
         (MatrixShape.Full vertB horizB heightB widthB) =
            MatrixShape.Full
               (ExtentPriv.Multiply vertA vertB)
               (ExtentPriv.Multiply horizA horizB)
               heightA widthB
   a <#> b =
      case unifyFactors (fullExtent a) (fullExtent b) of
         ((ExtentPriv.TagFact, ExtentPriv.TagFact), (unifyLeft, unifyRight)) ->
            multiply
               (mapExtent unifyLeft a)
               (mapExtent unifyRight b)

fullExtent ::
   Full vert horiz height width a ->
   Extent.Extent vert horiz height width
fullExtent = MatrixShape.fullExtent . Array.shape

unifyFactors ::
   (Extent.C vertA, Extent.C horizA, Extent.C vertB, Extent.C horizB) =>
   (ExtentPriv.Multiply vertA vertB ~ vertC) =>
   (ExtentPriv.Multiply horizA horizB ~ horizC) =>
   Extent.Extent vertA horizA height fuse ->
   Extent.Extent vertB horizB fuse width ->
   ((ExtentPriv.TagFact vertC, ExtentPriv.TagFact horizC),
    (Extent.Map vertA horizA vertC horizC height fuse,
     Extent.Map vertB horizB vertC horizC fuse width))
unifyFactors a b =
   ((ExtentPriv.multiplyTagLaw
         (ExtentPriv.heightFact a) (ExtentPriv.heightFact b),
     ExtentPriv.multiplyTagLaw
         (ExtentPriv.widthFact a) (ExtentPriv.widthFact b)),
    (ExtentPriv.Map $ flip ExtentPriv.unifyLeft b,
     ExtentPriv.Map $ ExtentPriv.unifyRight a))


instance
   (Extent.C vert, Extent.C horiz,
    Shape.C size, size ~ width, Eq width, Shape.C height) =>
      Multiply
         (MatrixShape.Full vert horiz height width)
         (MatrixShape.Hermitian size)
            where
   type Multiplied
         (MatrixShape.Full vert horiz height width)
         (MatrixShape.Hermitian size) =
            MatrixShape.Full vert horiz height width
   a <#> b = transpose $ Hermitian.multiplyFull Transposed b (transpose a)

instance
   (Extent.C vert, Extent.C horiz,
    Shape.C size, size ~ height, Eq height, Shape.C width) =>
      Multiply
         (MatrixShape.Hermitian size)
         (MatrixShape.Full vert horiz height width)
            where
   type Multiplied
         (MatrixShape.Hermitian size)
         (MatrixShape.Full vert horiz height width) =
            MatrixShape.Full vert horiz height width
   (<#>) = Hermitian.multiplyFull NonTransposed

instance
   (Shape.C shapeA, shapeA ~ shapeB, Eq shapeB) =>
      Multiply (MatrixShape.Hermitian shapeA) (MatrixShape.Hermitian shapeB)
         where
   type Multiplied
         (MatrixShape.Hermitian shapeA) (MatrixShape.Hermitian shapeB) =
            MatrixShape.Square shapeA
   a <#> b = Hermitian.multiplyFull NonTransposed a (Hermitian.toSquare b)


instance
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Extent.C vert, Extent.C horiz,
    Shape.C size, size ~ width, Eq width, Shape.C height) =>
      Multiply
         (MatrixShape.Full vert horiz height width)
         (MatrixShape.Triangular lo diag up size)
            where
   type Multiplied
         (MatrixShape.Full vert horiz height width)
         (MatrixShape.Triangular lo diag up size) =
            MatrixShape.Full vert horiz height width
   a <#> b =
      transpose $ Triangular.multiplyFull (Triangular.transpose b) (transpose a)

instance
   (MatrixShape.Content lo, MatrixShape.Content up, MatrixShape.TriDiag diag,
    Extent.C vert, Extent.C horiz,
    Shape.C size, size ~ height, Eq height, Shape.C width) =>
      Multiply
         (MatrixShape.Triangular lo diag up size)
         (MatrixShape.Full vert horiz height width)
            where
   type Multiplied
         (MatrixShape.Triangular lo diag up size)
         (MatrixShape.Full vert horiz height width) =
            MatrixShape.Full vert horiz height width
   (<#>) = Triangular.multiplyFull

instance
   (Shape.C sizeA, sizeA ~ sizeB, Eq sizeB,
    MultiplyTriangular loA upA loB upB,
    MatrixShape.TriDiag diagA, MatrixShape.TriDiag diagB) =>
      Multiply
         (MatrixShape.Triangular loA diagA upA sizeA)
         (MatrixShape.Triangular loB diagB upB sizeB) where
   type Multiplied
         (MatrixShape.Triangular loA diagA upA sizeA)
         (MatrixShape.Triangular loB diagB upB sizeB) =
            -- requires UndecidableInstances
            MultipliedTriangular loA diagA upA loB diagB upB sizeB
   (<#>) = multiplyTriangular

class
   (MatrixShape.Content loA, MatrixShape.Content upA,
    MatrixShape.Content loB, MatrixShape.Content upB) =>
      MultiplyTriangular loA upA loB upB where
   multiplyTriangular ::
      (Class.Floating a, Shape.C size, Eq size,
       MatrixShape.TriDiag diagA, MatrixShape.TriDiag diagB) =>
      Triangular loA diagA upA size a ->
      Triangular loB diagB upB size a ->
      Array (MultipliedTriangular loA diagA upA loB diagB upB size) a


type MultipliedTriangular loA diagA upA loB diagB upB size =
   ComposedTriangular
      (MultipliedPart loA loB)
      (MultipliedDiag diagA diagB)
      (MultipliedPart upA upB)
      size

type family MultipliedPart a b :: *
type instance MultipliedPart Empty b = b
type instance MultipliedPart Filled b = Filled

type family MultipliedDiag a b :: *
type instance MultipliedDiag Unit b = b
type instance MultipliedDiag NonUnit b = NonUnit

type family ComposedTriangular lo diag up size :: *
type instance ComposedTriangular Empty diag up size =
         MatrixShape.Triangular Empty diag up size
type instance ComposedTriangular Filled diag Empty size =
         MatrixShape.LowerTriangular diag size
type instance ComposedTriangular Filled diag Filled size =
         MatrixShape.Square size


instance MultiplyTriangular Empty Empty Empty Empty where
   multiplyTriangular = multiplyTriangularConform

instance MultiplyTriangular Empty Empty Filled Filled where
   multiplyTriangular a = Triangular.multiplyFull a . Triangular.toSquare

instance MultiplyTriangular Empty Filled Filled Filled where
   multiplyTriangular a = Triangular.multiplyFull a . Triangular.toSquare

instance MultiplyTriangular Filled Empty Filled Filled where
   multiplyTriangular a = Triangular.multiplyFull a . Triangular.toSquare

instance MultiplyTriangular Empty Filled Empty Filled where
   multiplyTriangular = multiplyTriangularConform

instance MultiplyTriangular Filled Empty Filled Empty where
   multiplyTriangular = multiplyTriangularConform

instance MultiplyTriangular Filled Empty Empty Filled where
   multiplyTriangular a = Triangular.multiplyFull a . Triangular.toSquare

instance MultiplyTriangular Empty Filled Filled Empty where
   multiplyTriangular a = Triangular.multiplyFull a . Triangular.toSquare

instance MultiplyTriangular Filled Filled Empty Empty where
   multiplyTriangular = multiplyTriangularToSquare

instance MultiplyTriangular Filled Filled Empty Filled where
   multiplyTriangular = multiplyTriangularToSquare

instance MultiplyTriangular Filled Filled Filled Empty where
   multiplyTriangular = multiplyTriangularToSquare

instance MultiplyTriangular Filled Filled Filled Filled where
   multiplyTriangular = multiplyTriangularToSquare

multiplyTriangularToSquare ::
   (MatrixShape.Content loA, MatrixShape.Content upA, MatrixShape.TriDiag diagA,
    MatrixShape.Content loB, MatrixShape.Content upB, MatrixShape.TriDiag diagB,
    Shape.C size, Eq size, Class.Floating a) =>
   Triangular loA diagA upA size a ->
   Triangular loB diagB upB size a ->
   Square size a
multiplyTriangularToSquare a b =
   transpose $ Triangular.multiplyFull (Triangular.transpose b) $
   transpose $ Triangular.toSquare a


newtype MultiplyTriangularConform lo up size a diagB diagA =
   MultiplyTriangularConform {
      getMultiplyTriangularConform ::
         Triangular lo diagA up size a ->
         Triangular lo diagB up size a ->
         Triangular lo (MultipliedDiag diagA diagB) up size a
   }

multiplyTriangularConform ::
   (Shape.C size, Eq size, Class.Floating a,
    MatrixShape.DiagUpLo lo up,
    MatrixShape.TriDiag diagA, MatrixShape.TriDiag diagB) =>
   (MultipliedDiag diagA diagB ~ diagC) =>
   Triangular lo diagA up size a ->
   Triangular lo diagB up size a ->
   Triangular lo diagC up size a
multiplyTriangularConform =
   getMultiplyTriangularConform $
   MatrixShape.switchTriDiag
      (MultiplyTriangularConform $ \a b ->
         Triangular.multiply (Triangular.relaxUnitDiagonal a) b)
      (MultiplyTriangularConform $ \a b ->
         Triangular.multiply a (Triangular.strictNonUnitDiagonal b))


instance
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vertA, Extent.C horizA,
    Extent.C vertB, Extent.C horizB,
    Shape.C heightA, Shape.C widthA, Shape.C widthB,
    widthA ~ heightB, Eq heightB) =>
      Multiply
         (MatrixShape.Full vertA horizA heightA widthA)
         (MatrixShape.Banded sub super vertB horizB heightB widthB)
            where
   type Multiplied
         (MatrixShape.Full vertA horizA heightA widthA)
         (MatrixShape.Banded sub super vertB horizB heightB widthB) =
            MatrixShape.Full
               (ExtentPriv.Multiply vertA vertB)
               (ExtentPriv.Multiply horizA horizB)
               heightA widthB
   a <#> b =
      case unifyFactors (fullExtent a) (bandedExtent b) of
         ((ExtentPriv.TagFact, ExtentPriv.TagFact), (unifyLeft, unifyRight)) ->
            transpose $
            Banded.multiplyFull
               (Banded.transpose $ Banded.mapExtent unifyRight b)
               (transpose $ mapExtent unifyLeft a)

instance
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vertA, Extent.C horizA,
    Extent.C vertB, Extent.C horizB,
    Shape.C heightA, Shape.C widthA, Shape.C widthB,
    widthA ~ heightB, Eq heightB) =>
      Multiply
         (MatrixShape.Banded sub super vertA horizA heightA widthA)
         (MatrixShape.Full vertB horizB heightB widthB)
            where
   type Multiplied
         (MatrixShape.Banded sub super vertA horizA heightA widthA)
         (MatrixShape.Full vertB horizB heightB widthB) =
            MatrixShape.Full
               (ExtentPriv.Multiply vertA vertB)
               (ExtentPriv.Multiply horizA horizB)
               heightA widthB
   a <#> b =
      case unifyFactors (bandedExtent a) (fullExtent b) of
         ((ExtentPriv.TagFact, ExtentPriv.TagFact), (unifyLeft, unifyRight)) ->
            Banded.multiplyFull
               (Banded.mapExtent unifyLeft a)
               (mapExtent unifyRight b)

instance
   (Unary.Natural subA, Unary.Natural superA,
    Unary.Natural subB, Unary.Natural superB,
    Extent.C vertA, Extent.C horizA,
    Extent.C vertB, Extent.C horizB,
    Shape.C heightA, Shape.C widthA, Shape.C widthB,
    widthA ~ heightB, Eq heightB) =>
      Multiply
         (MatrixShape.Banded subA superA vertA horizA heightA widthA)
         (MatrixShape.Banded subB superB vertB horizB heightB widthB) where
   type Multiplied
         (MatrixShape.Banded subA superA vertA horizA heightA widthA)
         (MatrixShape.Banded subB superB vertB horizB heightB widthB) =
            MatrixShape.Banded
               (subA :+: subB) (superA :+: superB)
               (ExtentPriv.Multiply vertA vertB)
               (ExtentPriv.Multiply horizA horizB)
               heightA widthB
   a <#> b =
      case unifyFactors (bandedExtent a) (bandedExtent b) of
         ((ExtentPriv.TagFact, ExtentPriv.TagFact), (unifyLeft, unifyRight)) ->
            Banded.multiply
               (Banded.mapExtent unifyLeft a)
               (Banded.mapExtent unifyRight b)

bandedExtent ::
   Banded.Banded sup super vert horiz height width a ->
   Extent.Extent vert horiz height width
bandedExtent = MatrixShape.bandedExtent . Array.shape


instance
   (Unary.Natural offDiag, Extent.C vert, Extent.C horiz,
    Shape.C size, size ~ width, Eq width, Shape.C height, Eq height) =>
      Multiply
         (MatrixShape.Full vert horiz height width)
         (MatrixShape.BandedHermitian offDiag size)
            where
   type Multiplied
         (MatrixShape.Full vert horiz height width)
         (MatrixShape.BandedHermitian offDiag size) =
            MatrixShape.Full vert horiz height width
   a <#> b = transpose $ BandedHermitian.multiplyFull Transposed b (transpose a)

instance
   (Unary.Natural offDiag, Extent.C vert, Extent.C horiz,
    Shape.C size, size ~ height, Eq height, Shape.C width, Eq width) =>
      Multiply
         (MatrixShape.BandedHermitian offDiag size)
         (MatrixShape.Full vert horiz height width)
            where
   type Multiplied
         (MatrixShape.BandedHermitian offDiag size)
         (MatrixShape.Full vert horiz height width) =
            MatrixShape.Full vert horiz height width
   (<#>) = BandedHermitian.multiplyFull NonTransposed

instance
   (Unary.Natural offDiag, Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz,
    Shape.C size, size ~ width, Eq width, Shape.C height, Eq height) =>
      Multiply
         (MatrixShape.Banded sub super vert horiz height width)
         (MatrixShape.BandedHermitian offDiag size)
            where
   type Multiplied
         (MatrixShape.Banded sub super vert horiz height width)
         (MatrixShape.BandedHermitian offDiag size) =
            MatrixShape.Banded
               (sub:+:offDiag) (super:+:offDiag) vert horiz height width
   a <#> b =
      Banded.multiply
         a (Banded.mapExtent Extent.fromSquare (BandedHermitian.toBanded b))

instance
   (Unary.Natural offDiag, Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz,
    Shape.C size, size ~ height, Eq height, Shape.C width, Eq width) =>
      Multiply
         (MatrixShape.BandedHermitian offDiag size)
         (MatrixShape.Banded sub super vert horiz height width)
            where
   type Multiplied
         (MatrixShape.BandedHermitian offDiag size)
         (MatrixShape.Banded sub super vert horiz height width) =
            MatrixShape.Banded
               (offDiag:+:sub) (offDiag:+:super) vert horiz height width
   a <#> b =
      Banded.multiply
         (Banded.mapExtent Extent.fromSquare (BandedHermitian.toBanded a)) b

instance
   (Unary.Natural offDiagA, Unary.Natural offDiagB,
    Shape.C sizeA, sizeA ~ sizeB, Shape.C sizeB, Eq sizeB) =>
      Multiply
         (MatrixShape.BandedHermitian offDiagA sizeA)
         (MatrixShape.BandedHermitian offDiagB sizeB)
            where
   type Multiplied
         (MatrixShape.BandedHermitian offDiagA sizeA)
         (MatrixShape.BandedHermitian offDiagB sizeB) =
            MatrixShape.Banded
               (offDiagA:+:offDiagB) (offDiagA:+:offDiagB)
               Small Small sizeA sizeB
   a <#> b =
      Banded.multiply (BandedHermitian.toBanded a) (BandedHermitian.toBanded b)
