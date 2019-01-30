{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE DeriveDataTypeable #-}

module Type.Data.Bool
    ( True
    , true
    , False
    , false
    , Not
    , not
    , (:&&:)
    , and
    , (:||:)
    , or
    , If
    , if_
    ) where

import Type.Base.Proxy (Proxy(Proxy))
import Data.Typeable (Typeable)

import qualified Prelude


data True deriving (Typeable)
true :: Proxy True
true = Proxy
instance Prelude.Show True where
    show _ = "True"
data False deriving (Typeable)
false :: Proxy False
false = Proxy
instance Prelude.Show False where
    show _ = "False"

type family Not x
type instance Not False = True
type instance Not True  = False
not :: Proxy x -> Proxy (Not x)
not Proxy = Proxy

type family x :&&: y
type instance False :&&: _x = False
type instance True  :&&: x = x
and :: Proxy x -> Proxy y -> Proxy (x :&&: y)
and Proxy Proxy = Proxy

type family x :||: y
type instance True  :||: _x = True
type instance False :||: x = x
or :: Proxy x -> Proxy y -> Proxy (x :||: y)
or Proxy Proxy = Proxy

type family If x y z
type instance If True y _z = y
type instance If False _y z = z
if_ :: Proxy x -> Proxy y -> Proxy z -> Proxy (If x y z)
if_ Proxy Proxy Proxy = Proxy
