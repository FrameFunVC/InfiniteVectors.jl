"""
    AbstractFiniteNZInfiniteVector{T} <: BiInfiniteVector{T}

Instance of BiInfiniteVector that implements `eachnonzeroindex(vec)`.

Also `hascompactsupport(vec) == true`
"""
abstract type AbstractFiniteNZInfiniteVector{T} <: BiInfiniteVector{T} end

"""
    hascompactsupport(vec::InfiniteVector)

Are their a finite number of non-zero elements.
"""
hascompactsupport(::AbstractFiniteNZInfiniteVector) = true
hascompactsupport(::InfiniteVector) = false

"""
    subvector(vec::CompactInfiniteVector) = vec.subvec

The vector of values at `eachnonzeroindex`
"""
@inline subvector(vec::AbstractFiniteNZInfiniteVector) = vec.subvec

"""
    sublength(vec::InfiniteVector) = vec.subvec

The number of non-zero indices
"""
@inline sublength(vec::AbstractFiniteNZInfiniteVector) = length(subvector(vec))

sum(vec::AbstractFiniteNZInfiniteVector) = sum(subvector(vec))

"""
    moment(vec::InfiniteVector, j)

The j-th discrete moment of a sequence is defined as \$\\sum_{k\\in ℤ} h_k k^j\$.
"""
function moment(vec::AbstractFiniteNZInfiniteVector, j)
    z = zero(eltype(vec))
    for k in eachnonzeroindex(vec)
        z += vec[k] * k^j
    end
    z
end

"""
    ztransform(vec::InfiniteVector, z)

The Z transform of a sequence is a continuous function of `z`, defined by
\$S(z) = \\sum_{k\\in ℤ} s_k z^{-k}\$
"""
function ztransform(vec::AbstractFiniteNZInfiniteVector, z)
    T = promote_type(eltype(vec), eltype(z), eltype(1/z))
    S = zero(T)
    for k in eachnonzeroindex(vec)
        S += vec[k] * z^T(-k)
    end
    S
end

for op in (:(==), :≈,)
    @eval $op(vec1::AbstractFiniteNZInfiniteVector, vec2::AbstractFiniteNZInfiniteVector) =
            $op(subvector(vec1), subvector(vec2)) && (eachnonzeroindex(vec1)==eachnonzeroindex(vec2))
end
