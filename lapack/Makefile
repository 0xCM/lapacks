run-test:
	runhaskell Setup configure --user --enable-tests -fbuildExamples
	runhaskell Setup build
	runhaskell Setup haddock
	./dist/build/lapack-test/lapack-test
	./dist/build/lapack-economic/lapack-economic
