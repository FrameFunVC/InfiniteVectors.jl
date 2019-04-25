
```@meta
DocTestSetup = quote
    using InfiniteVectors
end

```

# Indices
Even though an `InfiniteVector` is an `AbstractVector`
```jldoctest man
julia> InfiniteVector <: AbstractVector
true
```
it does not support indexing since
```jldoctest man
julia> vec = CompactInfiniteVector(1:3)
CompactInfiniteVector{Int64} with indices ℤ:
[  …, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, …  ]
julia> I = eachindex(vec)
ℤ
```
which is an object with infinite length.
```jldoctest man
julia> length(I)
∞
```

!!! note
    The function `eachnonzeroindex` is supported for vectors that have a finite number of non-zero elements.
    ```jldoctest man
    julia> eachnonzeroindex(vec)
    0:2
    ```
