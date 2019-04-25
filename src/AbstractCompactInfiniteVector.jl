"""
    AbstractFiniteNZInfiniteVector{T} <: AbstractDoubleInfiniteVector{T}

Instance of AbstractDoubleInfiniteVector that implements `eachnonzeroindex(vec)`.

Also `hascompactsupport(vec) == true`
"""
abstract type AbstractFiniteNZInfiniteVector{T} <: AbstractDoubleInfiniteVector{T} end

hascompactsupport(::AbstractFiniteNZInfiniteVector) = true
hascompactsupport(::AbstractDoubleInfiniteVector) = false

sum(vec::AbstractFiniteNZInfiniteVector) = sum(vec[i] for i in eachnonzeroindex)

function moment(vec::AbstractFiniteNZInfiniteVector, j)
    z = zero(eltype(vec))
    for k in eachnonzeroindex(vec)
        z += vec[k] * k^j
    end
    z
end

function ztransform(vec::AbstractFiniteNZInfiniteVector, z)
    T = promote_type(eltype(vec), eltype(z), eltype(1/z))
    S = zero(T)
    for k in eachnonzeroindex(vec)
        S += vec[k] * z^T(-k)
    end
    S
end
