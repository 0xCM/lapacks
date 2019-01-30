run:	run-blas run-lapack

run-blas:
	cabal run create-lapack-ffi -- \
	  data/blas.csv blas-src-3.7/ \
	  blas-ffi/src/Numeric/BLAS/FFI/ \
	  blas-comfort-array/src/Numeric/BLAS/ComfortArray/ \
	  blas-carray/src/Numeric/BLAS/CArray/

run-lapack:
	cabal run create-lapack-ffi -- \
	  --fortran-url='http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/%s' \
	  data/lapack.csv lapack-src-3.7/ \
	  lapack-ffi/src/Numeric/LAPACK/FFI/ \
	  lapack-comfort-array/src/Numeric/LAPACK/ComfortArray/ \
	  lapack-carray/src/Numeric/LAPACK/CArray/

# You may generate a lapack.list file from a listing of Fortran files.
# However, some of the functions seem to be unexported.
# E.g. there are several such functions with suffix 'x'.
# Instead you may extract the names from the library file:
#   nm /usr/lib/lapack/liblapack.a | fgrep T

csv:	data/blas-new.csv data/lapack-new.csv

data/%-new.csv:	data/%.list
	grep -v '^[ix]' $< | create-lapack-csv >$@
	grep '^[ix]' $< | create-lapack-csv >>$@
