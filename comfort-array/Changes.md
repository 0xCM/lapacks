# Change log for the `comfort-array` package

## 0.3

 * `Storable.Mutable.Array`: Replace `ForeignPtr` by `Array.Guarded.MutablePtr`.
   In the last release we altered the arrays after initialization
   which corrupted the debugging by the `guarded-allocation` package.
   This should be fixed now.

 * `Shape.sizeOffset`: It does not return a single offset anymore
   but an offset computation function.
   This allows to cache a size computation across many offset computations.

## 0.2

 * Add a monad parameter to the mutable `Storable` array type
   and generalize functions to `PrimMonad`s.
   This way the mutating functions can also be used in the `ST` monad.

## 0.1.2

 * Add immutable `Boxed` array type and mutable `Storable` array type.

## 0.1

 * Split `Shape.C` into `Shape.C` and `Shape.Indexed`.

## 0.0

 * Initial version featuring the `Shape.C` class with type function `Index`
   and the immutable `Storable` array type.
