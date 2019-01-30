test-build:
	runhaskell Setup configure --user --enable-tests
	runhaskell Setup build
	runhaskell Setup haddock

run-test:	test-build
	dist/build/comfort-array-test/comfort-array-test

%.html:	%.md
	pandoc --output=$@ $<
