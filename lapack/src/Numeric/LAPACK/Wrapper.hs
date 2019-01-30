module Numeric.LAPACK.Wrapper where

-- cf. Data.Bifunctor.Flip
newtype Flip f b a = Flip {getFlip :: f a b}
