abstract type AbstractFiniteNZInfiniteVector{T} <: AbstractDoubleInfiniteVector{T} end
"""
Indices of `vec` that are contain non-zero elements only possible if
```
hascompactsupport(vec) == true
```.
"""
function eachnonzeroindex(vec::AbstractDoubleInfiniteVector) end
hascompactsupport(::AbstractFiniteNZInfiniteVector) = true
hascompactsupport(::AbstractDoubleInfiniteVector) = false

sum(vec::AbstractFiniteNZInfiniteVector) = sum(vec[i] for i in eachnonzeroindex)

function moment(vec::AbstractFiniteNZInfiniteVector, j)
    z = zero(eltype(vec))
    for k in eachnonzeroindex(vec)
        z += s[k] * k^j
    end
    z
end

function ztransform(vec::AbstractFiniteNZInfiniteVector, z)
    T = promote_type(eltype(vec), eltype(z), eltype(1/z))
    S = zero(T)
    for k in eachnonzeroindex(s)
        S += s[k] * z^(-k)
    end
    S
end
