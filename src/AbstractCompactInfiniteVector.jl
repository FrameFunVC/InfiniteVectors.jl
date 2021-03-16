"""
    AbstractFiniteNZInfiniteVector{N,T} <: DoubleInfiniteArray{N,T} end

Instance of DoubleInfiniteArray that implements `eachnonzeroindex(vec)`.

Also `hascompactsupport(vec) == true`
"""
abstract type AbstractFiniteNZInfiniteArray{T,N} <: DoubleInfiniteArray{T,N} end
const AbstractFiniteNZInfiniteVector{T} = AbstractFiniteNZInfiniteArray{T,1}

"""
    hascompactsupport(vec::AbstractFiniteNZInfiniteArray)

Are their a finite number of non-zero elements.
"""
hascompactsupport(::AbstractFiniteNZInfiniteArray) = true
hascompactsupport(::DoubleInfiniteArray) = false
function Base.print_matrix(io::IO, X::AbstractFiniteNZInfiniteArray)
    Base.print_matrix(io,subarray(X))
end

"""
    subarray(A::AbstractFiniteNZInfiniteArray)

The vector of values at `eachnonzeroindex`
"""
@inline subarray(A::AbstractFiniteNZInfiniteArray) = A.subarray
@inline subvector(A::AbstractFiniteNZInfiniteVector) = subarray(A)

"""
    sublength(A::AbstractFiniteNZInfiniteArray)

The number of non-zero indices
"""
@inline sublength(A::AbstractFiniteNZInfiniteArray) = length(subarray(A))
"""
    subsize(A::AbstractFiniteNZInfiniteArray)

The number of non-zero indices
"""
@inline subsize(A::AbstractFiniteNZInfiniteArray) = size(subarray(A))

sum(A::AbstractFiniteNZInfiniteArray) = sum(subarray(A))

"""
    moment(vec::DoubleInfiniteVector, j)

The j-th discrete moment of a sequence is defined as \$\\sum_{k\\in ℤ} h_k k^j\$.
"""
function moment(vec::AbstractFiniteNZInfiniteVector, j)
    z = zero(eltype(vec))
    for k in eachnonzeroindex(vec)
        z += vec[k] * k.I[1]^j
    end
    z
end

"""
    ztransform(vec::DoubleInfiniteVector, z)

The Z transform of a sequence is a continuous function of `z`, defined by
\$S(z) = \\sum_{k\\in ℤ} s_k z^{-k}\$
"""
function ztransform(vec::AbstractFiniteNZInfiniteVector, z)
    T = promote_type(eltype(vec), eltype(z), eltype(1/z))
    S = zero(T)
    for k in eachnonzeroindex(vec)
        S += vec[k] * z^T(-k.I[1])
    end
    S
end

for op in (:(==), :≈,)
    @eval $op(vec1::AbstractFiniteNZInfiniteArray, vec2::AbstractFiniteNZInfiniteArray) =
            $op(subarray(vec1), subarray(vec2)) && (eachnonzeroindex(vec1)==eachnonzeroindex(vec2))
end

Base.axes(A::AbstractFiniteNZInfiniteArray{T,N}) where {T,N} = ntuple(k->subsize(A)[k]==1 ? (0:0) : Integers(),Val(N))
Base.size(A::AbstractFiniteNZInfiniteArray{T,N}) where {T,N} = ntuple(k->subsize(A)[k]==1 ? 1 : DoubleInfinity(),Val(N))
