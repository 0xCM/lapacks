{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}
module Numeric.LAPACK.Format (
   (##),
   Format(format),
   FormatArray(formatArray),
   deflt,
   ) where

import qualified Numeric.LAPACK.Matrix.Shape.Private as MatrixShape
import qualified Numeric.LAPACK.Matrix.Extent.Private as Extent
import Numeric.LAPACK.Matrix.Shape.Private
         (Order(RowMajor, ColumnMajor), Filled(Filled), UnaryProxy)
import Numeric.LAPACK.Matrix.Private (Full)
import Numeric.LAPACK.Scalar (conjugate)
import Numeric.LAPACK.Wrapper (Flip(Flip, getFlip))

import qualified Numeric.Netlib.Class as Class

import qualified Type.Data.Num.Unary.Literal as TypeNum
import qualified Type.Data.Num.Unary as Unary
import Type.Data.Num (integralFromProxy)

import qualified Data.Array.Comfort.Storable.Unchecked as Array
import qualified Data.Array.Comfort.Shape as Shape
import Data.Array.Comfort.Storable (Array)

import qualified Text.PrettyPrint.Boxes as TextBox
import Text.PrettyPrint.Boxes (Box, (/+/))
import Text.Printf (PrintfArg, printf)

import qualified Data.List.Reverse.StrictSpine as ListRev
import qualified Data.List.Match as Match
import qualified Data.List.HT as ListHT
import qualified Data.List as List
import Data.Functor.Compose (Compose(Compose, getCompose))
import Data.Foldable (foldMap)
import Data.List (mapAccumL, transpose)
import Data.Complex (Complex((:+)))
import Data.Maybe.HT (toMaybe)
import Data.Maybe (fromMaybe)
import Data.Char (isSpace)


infix 0 ##

(##) :: (Format a) => a -> String -> IO ()
a ## fmt = putStr $ trim $ TextBox.render $ format fmt a

trim :: String -> String
trim = unlines . map (ListRev.dropWhile isSpace) . lines


deflt :: String
deflt = "%.4g"


class Format a where
   format :: String -> a -> Box

instance Format Int where
   format _fmt = TextBox.text . show

instance Format Float where
   format fmt = TextBox.text . printf fmt

instance Format Double where
   format fmt = TextBox.text . printf fmt

instance (Class.Real a) => Format (Complex a) where
   format fmt = TextBox.text . concat . printfComplex fmt

instance (Format a) => Format [a] where
   format fmt = TextBox.vsep 1 TextBox.right . map (format fmt)

instance (Format a, Format b) => Format (a,b) where
   format fmt (a,b) = format fmt a /+/ format fmt b

instance (Format a, Format b, Format c) => Format (a,b,c) where
   format fmt (a,b,c) = format fmt a /+/ format fmt b /+/ format fmt c

instance (FormatArray sh, Class.Floating a) => Format (Array sh a) where
   format = formatArray


class (Shape.C sh) => FormatArray sh where
   {-
   We use constraint @(Class.Floating a)@ and not @(Format a)@
   because it allows us to align the components of complex numbers.
   -}
   formatArray :: (Class.Floating a) => String -> Array sh a -> Box

instance (Integral i) => FormatArray (Shape.ZeroBased i) where
   formatArray = formatVector

instance (Integral i) => FormatArray (Shape.OneBased i) where
   formatArray = formatVector

formatVector :: (Shape.C sh, Class.Floating a) => String -> Array sh a -> Box
formatVector fmt =
   TextBox.hsep 1 TextBox.right .
   map (TextBox.text . concat . printfFloating fmt) . Array.toList

instance
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
      FormatArray (MatrixShape.Full vert horiz height width) where
   formatArray = formatFull

formatFull ::
   (Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   String -> Full vert horiz height width a -> Box
formatFull fmt m =
   let MatrixShape.Full order extent = Array.shape m
   in  formatAligned (printfFloating fmt) $
       splitRows order (Extent.dimensions extent) $ Array.toList m

instance
   (Eq lower, Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
      FormatArray (MatrixShape.Split lower vert horiz height width) where
   formatArray = formatHouseholder

formatHouseholder ::
   (Eq lower, Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width,
    Class.Floating a) =>
   String -> Array (MatrixShape.Split lower vert horiz height width) a -> Box
formatHouseholder fmt m =
   let MatrixShape.Split _ order extent = Array.shape m
   in formatSeparateTriangle (printfFloating fmt) $
      splitRows order (Extent.dimensions extent) $ Array.toList m

instance (Shape.C size) => FormatArray (MatrixShape.Hermitian size) where
   formatArray = formatHermitian

formatHermitian ::
   (Shape.C size, Class.Floating a) =>
   String -> Array (MatrixShape.Hermitian size) a -> Box
formatHermitian fmt m =
   let MatrixShape.Hermitian order size = Array.shape m
   in  formatSeparateTriangle (printfFloating fmt) $
       complementTriangle conjugate order (Shape.size size) $ Array.toList m

formatSymmetric ::
   (Shape.C size, Class.Floating a) =>
   String -> Array (MatrixShape.Symmetric size) a -> Box
formatSymmetric fmt m =
   let MatrixShape.Triangular _diag (Filled, Filled) order size = Array.shape m
   in  formatSeparateTriangle (printfFloating fmt) $
       complementTriangle id order (Shape.size size) $ Array.toList m

complementTriangle ::
   (Class.Floating a) => (a -> a) -> Order -> Int -> [a] -> [[a]]
complementTriangle adapt order n xs =
   let mergeTriangles lower upper =
         zipWith (++) (map (map adapt . init) lower) upper
   in case order of
         RowMajor ->
            let tri = slice (take n $ iterate pred n) xs
                trans = reverse $ transpose $ map reverse tri
            in  mergeTriangles trans tri
         ColumnMajor ->
            let tri = slice (take n [1..]) xs
            in  mergeTriangles tri (transpose tri)

instance
   (MatrixShape.Content lo, MatrixShape.Content up,
    MatrixShape.TriDiag diag, Shape.C size) =>
      FormatArray (MatrixShape.Triangular lo diag up size) where
   formatArray fmt =
      getFormatTriangular $
      MatrixShape.switchDiagUpLoSym
         (FormatTriangular $ \m ->
            let MatrixShape.Triangular _diag _uplo order size = Array.shape m
                n0 = Unary.unary TypeNum.u0
            in formatAligned (printfFloatingMaybe fmt) $
               formatBanded (n0,n0) order (size,size) $ Array.toList m)
         (FormatTriangular $ formatTriangular fmt)
         (FormatTriangular $ formatTriangular fmt)
         (FormatTriangular $
            formatSymmetric fmt .
            Array.mapShape MatrixShape.strictNonUnitDiagonal)

newtype FormatTriangular diag size a b lo up =
   FormatTriangular {
      getFormatTriangular ::
         Array (MatrixShape.Triangular lo diag up size) a -> b
   }

formatTriangular ::
   (MatrixShape.TriDiag diag, MatrixShape.UpLo lo up,
    Shape.C size, Class.Floating a) =>
   String -> Array (MatrixShape.Triangular lo diag up size) a -> Box
formatTriangular fmt m =
   let MatrixShape.Triangular _diag uplo order size = Array.shape m
   in  formatAligned (printfFloatingMaybe fmt) $
       MatrixShape.caseLoUp uplo
         padLowerTriangle padUpperTriangle order (Shape.size size) $
       Array.toList m

padUpperTriangle :: Order -> Int -> [a] -> [[Maybe a]]
padUpperTriangle order n xs =
   let mxs = map Just xs
       nothings = iterate (Nothing:) []
   in case order of
         RowMajor ->
            zipWith (++) nothings (slice (take n $ iterate pred n) mxs)
         ColumnMajor ->
            transpose $
            zipWith (++)
               (slice (take n [1..]) mxs)
               (reverse $ take n nothings)

padLowerTriangle :: Order -> Int -> [a] -> [[Maybe a]]
padLowerTriangle order n xs =
   map (map Just) $
   case order of
      RowMajor -> slice (take n [1..]) xs
      ColumnMajor ->
         foldr (\(y:ys) zs -> [y] : zipWith (:) ys zs) []
            (slice (take n $ iterate pred n) xs)

slice :: [Int] -> [a] -> [[a]]
slice ns xs =
   snd $ mapAccumL (\ys n -> let (vs,ws) = splitAt n ys in (ws,vs)) xs ns

formatSeparateTriangle :: (a -> [String]) -> [[a]] -> Box
formatSeparateTriangle printFmt =
   alignSeparated . map concat .
   zipWith
      (zipWith (\sep -> attachSeparators sep . printFmt))
      (ListHT.outerProduct
         (\row col -> if row==col then Bar else Space)
         [(0::Int)..] [0..])


instance
   (Unary.Natural sub, Unary.Natural super,
    Extent.C vert, Extent.C horiz, Shape.C height, Shape.C width) =>
      FormatArray (MatrixShape.Banded sub super vert horiz height width) where
   formatArray fmt m =
      let MatrixShape.Banded offDiag order extent = Array.shape m
      in  formatAligned (printfFloatingMaybe fmt) $
          formatBanded offDiag order (Extent.dimensions extent) $
          Array.toList m

formatBanded ::
   (Shape.C height, Shape.C width, Unary.Natural sub, Unary.Natural super) =>
   (UnaryProxy sub, UnaryProxy super) -> Order ->
   (height, width) -> [a] -> [[Maybe a]]
formatBanded (sub,super) order (height,width) xs =
   let slices =
         ListHT.sliceVertical (MatrixShape.bandedBreadth (sub,super)) xs
       m = Shape.size height
       n = Shape.size width
   in case order of
         RowMajor ->
            map (take n) $
            zipWith (shiftRow Nothing)
               (iterate (1+) (- integralFromProxy sub))
               (map (map Just) slices)
         ColumnMajor ->
            let ku = integralFromProxy super
            in take m $ drop ku $
               foldr
                  (\col mat ->
                     zipWith (:) (map Just col ++ repeat Nothing) ([]:mat))
                  (replicate (ku + m - n) [])
                  slices


instance
   (Unary.Natural offDiag, Shape.C size) =>
      FormatArray (MatrixShape.BandedHermitian offDiag size) where
   formatArray fmt m =
      let MatrixShape.BandedHermitian offDiag order size = Array.shape m
      in  formatSeparateTriangle (printfFloatingMaybe fmt) $
          formatBandedHermitian offDiag order size $ Array.toList m

formatBandedHermitian ::
   (Unary.Natural offDiag, Shape.C size, Class.Floating a) =>
   UnaryProxy offDiag -> Order -> size -> [a] -> [[Maybe a]]
formatBandedHermitian offDiag order _size xs =
   let k = integralFromProxy offDiag
       slices = ListHT.sliceVertical (k + 1) xs
   in case order of
         RowMajor ->
            foldr
               (\row square ->
                  Match.take ([]:square) (map Just row)
                  :
                  zipWith (:)
                     (tail $ map (Just . conjugate) row ++ repeat Nothing)
                     square)
               [] slices
         ColumnMajor ->
            zipWith (shiftRow Nothing) (iterate (1+) (-k)) $ map (map Just) $
            zipWith (++)
               (map (map conjugate . init) slices)
               (drop k $
                foldr
                  (\column band ->
                     zipWith (++) (map (:[]) column ++ repeat []) ([]:band))
                  (replicate k [])
                  slices)

shiftRow :: a -> Int -> [a] -> [a]
shiftRow pad k = if k<=0 then drop (-k) else (replicate k pad ++)

splitRows ::
   (Shape.C height, Shape.C width) =>
   Order -> (height, width) -> [a] -> [[a]]
splitRows order (height,width) =
   case order of
      RowMajor -> ListHT.sliceVertical (Shape.size width)
      ColumnMajor -> ListHT.sliceHorizontal (Shape.size height)

formatAligned :: (a -> [String]) -> [[a]] -> Box
formatAligned printFmt =
   alignSeparated . map (concatMap (attachSeparators Space . printFmt))


data Separator = Empty | Space | Bar
   deriving (Eq, Ord, Show)

alignSeparated :: [[(Separator, String)]] -> Box
alignSeparated =
   TextBox.hcat TextBox.top .
   map (TextBox.vcat TextBox.right . map TextBox.text) .
   concatMap ((\(seps,column) -> [map formatSeparator seps, column]) . unzip) .
   List.unfoldr (viewLAll (Empty,""))

viewLAll :: a -> [[a]] -> Maybe ([a], [[a]])
viewLAll x0 xs =
   toMaybe (any (not.null) xs)
      (unzip $ map (fromMaybe (x0,[]) . ListHT.viewL) xs)

formatSeparator :: Separator -> String
formatSeparator sep = case sep of Empty -> ""; Space -> " "; Bar -> "|"

attachSeparators :: Separator -> [str] -> [(Separator, str)]
attachSeparators sep = zip (sep:repeat Empty)


printfFloating :: (Class.Floating a) => String -> a -> [String]
printfFloating fmt =
   getFlip $
   Class.switchFloating
      (Flip $ (:[]) . printf fmt)
      (Flip $ (:[]) . printf fmt)
      (Flip $ printfComplex fmt)
      (Flip $ printfComplex fmt)

printfFloatingMaybe :: (Class.Floating a) => String -> Maybe a -> [String]
printfFloatingMaybe fmt =
   getFlip $ getCompose $
   Class.switchFloating
      (Compose $ Flip $ (:[]) . foldMap (printf fmt))
      (Compose $ Flip $ (:[]) . foldMap (printf fmt))
      (Compose $ Flip $ maybe ["",""] (printfComplex fmt))
      (Compose $ Flip $ maybe ["",""] (printfComplex fmt))

printfComplex :: (Class.Real a) => String -> Complex a -> [String]
printfComplex fmt =
   getFlip $ getCompose $
   Class.switchReal
      (Compose $ Flip $ printfComplexAux fmt)
      (Compose $ Flip $ printfComplexAux fmt)

printfComplexAux ::
   (PrintfArg a, Class.Real a) => String -> Complex a -> [String]
printfComplexAux fmt (r:+i) =
   if i<0 || isNegativeZero i
     then [printf (fmt ++ "-") r, printf (fmt ++ "i") (-i)]
     else [printf (fmt ++ "+") r, printf (fmt ++ "i") i]
