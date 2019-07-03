"""
    mutable struct CompactInfiniteVector{T} <: AbstractFiniteNZInfiniteVector{T}

A CompactInfiniteVector contains a sequence of nonzero elements starting at a given offset.

    CompactInfiniteVector(a::AbstractVector{T}, offset = 0)

Construct a CompactInfiniteVector with indices starting at offset.

# Examples
```jldoctest; setup = :(using InfiniteVectors)
julia> CompactInfiniteVector(1:3,-1)
CompactInfiniteVector{Int64} with indices ℤ:
[  …, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, …  ]
```
"""
mutable struct CompactInfiniteVector{T} <: AbstractFiniteNZInfiniteVector{T}
    subvec  ::  Vector{T}
    offset  ::  Int

    CompactInfiniteVector{T}(a, offset) where T = new(a, offset)
end

"""
    δ(n::Int)

Construct Dirac delts, i.e., a bi-infinite vector with δ(k)=1, if k=n, and δ(k)=0, otherwise.
"""
δ(k::Int) = CompactInfiniteVector([1], k)


CompactInfiniteVector(a::Vector{T}, offset = 0) where {T} = CompactInfiniteVector{T}(a, offset)
CompactInfiniteVector(a::AbstractVector{T}, offset = 0) where {T} = CompactInfiniteVector{T}(collect(a), offset)

@inline offset(vec::CompactInfiniteVector) = vec.offset



@inline setindex!(vec::CompactInfiniteVector, k::Int) =
(k < offset(vec)) || (k >= offset(vec)+sublength(vec)) ? error("Not possible") : setindex!(subvector(vec), _mapindex(vec, k))


"""
    shift!(vec::InfiniteArrays, m::Int)

In-place shifting of vector. (not always possible)
"""
shift!(vector::CompactInfiniteVector, k::Int) = (vector.offset = offset(vector)+k; vector)


"""
    reverse!(vec::CompactInfiniteVector)

In-place time reversel: `vec(-k)`
"""
reverse!(vec::CompactInfiniteVector) = (reverse!(subvector(vec)); vec.offset = -_lastindex(vec);vec)



"""
    struct FixedInfiniteVector{L,OFS,T} <: AbstractFiniteNZInfiniteVector{T}

A FixedInfiniteVector is a fully typed sequence of nonzero elements starting at a given offset.

    FixedInfiniteVector(a::AbstractVector{T}, offset = 0)

Construct a FixedInfiniteVector with indices starting at offset.

# Examples
```jldoctest; setup = :(using InfiniteVectors)
julia> FixedInfiniteVector(1:3,-1)
FixedInfiniteVector{3,-1,Int64} with indices ℤ:
[  …, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, …  ]
```
"""
struct FixedInfiniteVector{L,OFS,T} <: AbstractFiniteNZInfiniteVector{T}
    subvec       ::  SVector{L,T}

    FixedInfiniteVector(a::SVector{L,T}, ::Val{OFS}) where {L,T,OFS} = new{L,OFS,T}(a)
end


# type unstable constructor
"""
    FixedInfiniteVector(a::AbstractVector{T}, offset = 0)

Construct a FixedInfiniteVector with indices starting at offset.
"""
FixedInfiniteVector(a::AbstractVector, offset = 0) = FixedInfiniteVector(length(a)==0 ? SVector{0,eltype(a)}() : SVector(a...), Val{offset}())

# type unstable constructor
FixedInfiniteVector(vec::CompactInfiniteVector) = FixedInfiniteVector(subvector(vec), offset(vec))

CompactInfiniteVector(vec::FixedInfiniteVector{L,OFS,T}) where {L,OFS,T} = CompactInfiniteVector(Array(subvector(vec)), OFS)

offset(vec::FixedInfiniteVector{L,OFS}) where {L,OFS} = OFS

sublength(vec::FixedInfiniteVector{L,OFS,T}) where {L,OFS,T} = L






"The first even number greater than or equal to n."
nexteven(n) = isodd(n) ? n+1 : n

"The last even number, smaller than or equal to n."
previouseven(n) = isodd(n) ? n-1 : n

"The first odd number greater than or equal to n."
nextodd(n) = isodd(n) ? n : n+1





for COMPACTVECTOR in (:CompactInfiniteVector,:FixedInfiniteVector)
    @eval function copy(bc::Base.Broadcast.Broadcasted{<:Base.Broadcast.ArrayStyle{<:$COMPACTVECTOR},<:Tuple,F,Tuple{T}} where {F,T<:$COMPACTVECTOR})
        vc = bc.args[1]
        CompactInfiniteVector((bc.f).(subvector(vc)), offset(vc))
    end

    @eval copy(vec::$COMPACTVECTOR) = $COMPACTVECTOR(copy(subvector(vec)), offset(vec))

    @eval eachnonzeroindex(vec::$COMPACTVECTOR) = (1:sublength(vec)) .+ (offset(vec)-1)

    @eval shift(vector::$COMPACTVECTOR, k::Int) = $COMPACTVECTOR(copy(subvector(vector)), offset(vector)+k)

    @eval reverse(vec::$COMPACTVECTOR) = $COMPACTVECTOR(reverse(subvector(vec), 1), -_lastindex(vec))

    @eval downsample(vec::$COMPACTVECTOR, m::Int) =
        $COMPACTVECTOR(subvector(vec)[mod(m-offset(vec), m)+1:m:end], cld(offset(vec), m))

    @eval function upsample(vec::$COMPACTVECTOR, m::Int)
        v = similar(vec, m*(sublength(vec)-1)+1)
        fill!(v, 0)
        v[1:m:end] .= subvector(vec)
        $COMPACTVECTOR(v, m*offset(vec))
    end

    @eval @inline _mapindex(vec::$COMPACTVECTOR, k) = k - offset(vec) + 1
    @eval @inline i_mapindex(vec::$COMPACTVECTOR, l) = l + offset(vec) - 1
    @eval @inline _firstindex(vec::$COMPACTVECTOR) = i_mapindex(vec, 1)
    @eval @inline _lastindex(vec::$COMPACTVECTOR) = i_mapindex(vec, sublength(vec))

    @eval @inline getindex(vec::$COMPACTVECTOR, k::Int) =
        (k < offset(vec)) || (k >= offset(vec)+sublength(vec)) ? zero(eltype(vec)) : getindex(subvector(vec), _mapindex(vec, k))

    @eval function getindex(vec::$COMPACTVECTOR, c::UnitRange)
        e = zeros(eltype(vec), length(c))
        i1 = max(offset(vec), c[1])
        i2 = min(offset(vec) + sublength(vec)-1, c[end])
        e1 = max(1, 1-c[1]+offset(vec))
        # e2 = min(length(e),i2-i1+1-c[1]+offset(vec))
        e[(0:(i2-i1)) .+ e1] .= subvector(vec)[(i1:i2) .+ (-offset(vec)+1)]
        e
    end

    @eval function alternating_flip(vec::$COMPACTVECTOR, pivot = 1)
        hflip = similar(subvector(vec))
        # The first element of hflip has index pivot-_lastindex(f).
        # Whether or not we need to flip its sign depend on the parity of this index:
        isodd(pivot-_lastindex(vec)) ? t = -1 : t = 1
        hflip[1:2:end] =  t * subvector(vec)[end:-2:1]
        hflip[2:2:end] = -t * subvector(vec)[end-1:-2:1]
        $COMPACTVECTOR(hflip, pivot - _lastindex(vec))
    end

    @eval function alternating(vec::$COMPACTVECTOR)
        halt = similar(subvector(vec))
        t = (-1)^_firstindex(vec)
        for k in eachindex(subvector(vec))
            halt[k] = t * subvector(vec)[k]
            t = -t
        end
        $COMPACTVECTOR(halt, _firstindex(vec))
    end

    @eval evenpart(vec::$COMPACTVECTOR) =
        $COMPACTVECTOR([vec[j] for j in nexteven(_firstindex(vec)):2:_lastindex(vec)], div(nexteven(_firstindex(vec)),2))

    @eval oddpart(vec::$COMPACTVECTOR) =
        $COMPACTVECTOR([vec[j] for j in nextodd(_firstindex(vec)):2:_lastindex(vec)], div(previouseven(_firstindex(vec)),2))

    @eval support(vec::$COMPACTVECTOR) = (offset(vec), offset(vec) + sublength(vec) - 1)

    @eval support(vec::$COMPACTVECTOR, j::Int, k::Int) = (1/(1<<j)*(support(vec)[1]+k), 1/(1<<j)*(support(vec)[2]+k))

    @eval function inv(a::$COMPACTVECTOR, m::Int; R=sublength(a))
        T = eltype(a)
        if m == 1
            return inv(a)
        end
        n = sublength(a)
        iseven(n) && (n += 1)
        iseven(R) ? R == R>>1 : (R-1) >>1
        l = n >> 1
        R = max(l,R)
        # Symmetrise such that I ≈ -n/2..n/2
        sym_shift = -l-offset(a)
        b = shift(a, sym_shift)
        I = -R:R
        r = cld(2R,m)

        # there are some rows, therefore, limit Ii
        Ii = -r:r
        Ij = I
        A = zeros(T, length(Ii), length(Ij))
        for i in Ii
            for j in Ij
                A[i-Ii[1]+1, j-Ij[1]+1] = b[m*i-j]
            end
        end
        # determine rhs
        e = δ(0)[Ii]

        # solve and shift back
        $COMPACTVECTOR(pinv(A)*e, sym_shift+first(I))
    end

end

@eval conv(vec1::Union{CompactInfiniteVector,FixedInfiniteVector}, vec2::Union{CompactInfiniteVector,FixedInfiniteVector}) =
    CompactInfiniteVector(conv(convert(Vector, subvector(vec1)), convert(Vector, subvector(vec2))), offset(vec1) + offset(vec2))

@eval conv(vec1::FixedInfiniteVector, vec2::FixedInfiniteVector) =
    FixedInfiniteVector(conv(convert(Vector,subvector(vec1)), convert(Vector,subvector(vec2))), offset(vec1) + offset(vec2))



@doc """
    reverse(vec::BiInfiniteVector)

Time-reversel: `vec(-k)`
"""
reverse(vec)

@doc """
    alternating_flip(vec::BiInfiniteVector, pivot = 1)

From a given filter \$h(i)\$, compute a new filter satisfying the alternating flip relation,
centered around the given pivot:

``g(k) = (-1)^k h(pivot-k)``
"""
alternating_flip(vec, pivot = 1)

@doc """
    alternating(vec::InfiniteVector)

Compute a new 'alternating' filter satisfying \$g(k) = (-1)^k h(k)`.\$
"""
alternating(vec)

@doc """
evenpart(vec::InfiniteVector)

    Return the even part of a sequence `s`, defined by `s_e[k] = s[2k]`.
"""
evenpart(vec)

@doc """
    oddpart(vec::InfiniteVector)

Return the odd part of a sequence `s`, defined by `s_o[k] = s[2k+1]`
"""
oddpart(vec)

@doc """
    inv(a::CompactInfiniteVector{T}, m::Int)

A solution of \$[a*b]_{↓m}=δ\$, given a.
"""
inv(vec::CompactInfiniteVector)

@doc """
    eachnonzeroindex(vec)

Indices of `vec` that contain non-zero elements. This function can be called if

```
hascompactsupport(vec) == true
```
"""
eachnonzeroindex(vec)

@doc """
    support(vec::CompactInfiniteVector)

The minimum and maximum index of the non-zero elements
"""
support(vec)
