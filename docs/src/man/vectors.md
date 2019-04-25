# Vector types

This sections contains a list of the different supported vector types accompanied with the functions specific to those types.

```@contents
Pages = ["vectors.md"]
Depth = 2
```

## Compact vectors

The abstract type for vectors with a finite number of non-zero elements is

```@docs
InfiniteVectors.AbstractFiniteNZInfiniteVector
```

### CompactInfiniteVector
There are two types that support compact vectors. The second one is fully typed.
```@docs
CompactInfiniteVector
FixedInfiniteVector
```
For the dirac delta a special constructor is present. 
```@docs
δ
```

### Functions specific to compact vectors
```@docs
eachnonzeroindex
subvector
support
```

## Periodic vectors
```@docs
InfiniteVectors.AbstractPeriodicInfiniteVector
```
```@docs
PeriodicInfiniteVector
```

### Functions specific to periodic vectors
```@docs
period
```
