test:
	runhaskell Setup configure --user --enable-tests --enable-benchmarks
	runhaskell Setup build
	runhaskell Setup haddock
	./dist/build/test/test
	(cd tests; sh runtests.sh)
