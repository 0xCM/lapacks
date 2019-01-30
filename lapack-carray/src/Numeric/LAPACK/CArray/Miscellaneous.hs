-- Do not edit! Automatically generated by create-lapack-ffi.
module Numeric.LAPACK.CArray.Miscellaneous where

import qualified Numeric.LAPACK.FFI.Miscellaneous as FFI
import qualified Numeric.Netlib.CArray.Utility as Call

import Foreign.C.Types (CInt, CChar)

import Control.Monad.Trans.Cont (evalContT)
import Control.Monad.IO.Class (liftIO)


-- | <http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chla_transtype.f>
chla_transtype ::
   Int {- ^ trans -} ->
   IO CChar
chla_transtype trans = do
   evalContT $ do
      transPtr <- Call.cint trans
      liftIO $ FFI.chla_transtype transPtr

-- | <http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f>
ieeeck ::
   Int {- ^ ispec -} ->
   Float {- ^ zero -} ->
   Float {- ^ one -} ->
   IO CInt
ieeeck ispec zero one = do
   evalContT $ do
      ispecPtr <- Call.cint ispec
      zeroPtr <- Call.float zero
      onePtr <- Call.float one
      liftIO $ FFI.ieeeck ispecPtr zeroPtr onePtr

-- | <http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladiag.f>
ladiag ::
   Char {- ^ diag -} ->
   IO CInt
ladiag diag = do
   evalContT $ do
      diagPtr <- Call.char diag
      liftIO $ FFI.ladiag diagPtr

-- | <http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaprec.f>
laprec ::
   Char {- ^ prec -} ->
   IO CInt
laprec prec = do
   evalContT $ do
      precPtr <- Call.char prec
      liftIO $ FFI.laprec precPtr

-- | <http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilatrans.f>
latrans ::
   Char {- ^ trans -} ->
   IO CInt
latrans trans = do
   evalContT $ do
      transPtr <- Call.char trans
      liftIO $ FFI.latrans transPtr

-- | <http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilauplo.f>
lauplo ::
   Char {- ^ uplo -} ->
   IO CInt
lauplo uplo = do
   evalContT $ do
      uploPtr <- Call.char uplo
      liftIO $ FFI.lauplo uploPtr
