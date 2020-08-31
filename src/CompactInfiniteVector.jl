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
mutable struct CompactInfiniteArray{T,N,A<:Array{T,N}} <: AbstractFiniteNZInfiniteArray{T,N}
    subarray  ::  A
    offset  ::  NTuple{N,Int}

    CompactInfiniteArray{T,N}(a::AbstractArray{T,N}, offset::NTuple{N,Int}) where {T,N} = new{T,N,typeof(a)}(a, offset)
    CompactInfiniteArray{T,1}(a::AbstractArray{T,N}, offset::Int) where {T,N} = CompactInfiniteArray{T,1}(a, tuple(offset))
    CompactInfiniteArray(a, offset) where N = CompactInfiniteArray{eltype(a),ndims(a)}(a, offset)
end
const CompactInfiniteVector{T} = CompactInfiniteArray{T,1}

function Broadcast.materialize(bc::Broadcast.Broadcasted{<:InfiniteArrayStyle{<:AbstractFiniteNZInfiniteArray}})
    ixs = eachnonzeroindex.(bc.args)
    i1s = _bcfirst(bc.args...)
    ins = _bclast(bc.args...)
    sarrays = map(x->_resize(ins.-i1s.+1,i1s,__first(x),subarray(x)),bc.args)
    CompactInfiniteArray(Broadcast.broadcast(bc.f,  sarrays...), i1s)
end
_equal(s1,s2) = ((s1[1]==1) || (s1[1]==s2[1])) && _equal(Base.tail(s1),Base.tail(s2))
_equal(::Tuple{},s2) = true

function _resize(s, first, _first, array)
    if _equal(size(array),s)
        return array
    end
    A = zeros(eltype(array),s)
    I = eachindex(IndexCartesian(),array)
    copyto!(A,I.+CartesianIndex(_fill(first,_first)),array,I)
    A
end
_fill(t1::NTuple{N,Integer},t2::Tuple{Vararg{Union{Integer,Nothing}}}) where N =
    ntuple(k->isnothing(t2[k]) ? 0 : t2[k]-t1[k] ,Val(length(t2)))
function __first(A::AbstractFiniteNZInfiniteArray)
    os = _offset(A)
    ndims(A)==1 && (return os)
    ntuple(k->size(A,k)==1 ? nothing : os[k] ,Val(ndims(A)))
end
function __last(A::AbstractFiniteNZInfiniteArray)
    os = Tuple(last(eachnonzeroindex(A)))
    ndims(A)==1 && (return os)
    ntuple(k->size(A,k)==1 ? nothing : os[k] ,Val(ndims(A)))
end
for (_bcfirst,__first,broadcast_first,_bcs,_bcs1,_m) in (
        (:_bcfirst,:__first,:broadcast_first,:_bcs, :_bcs1, :min),
        (:_bclast,  :__last,:broadcast_last, :_bcsl,:_bcs1l,:max),
    )
    @eval begin
        $_bcfirst(A, B...) = $broadcast_first($__first(A), $_bcfirst(B...))
        $_bcfirst(A) = $__first(A)
        $broadcast_first(shape) = shape
        $broadcast_first(shape, shape1, shapes...) = $broadcast_first($_bcs(shape, shape1), shapes...)
        $_bcs(::Tuple{}, ::Tuple{}) = ()
        $_bcs(::Tuple{}, newshape::Tuple) = (newshape[1], $_bcs((), Base.tail(newshape))...)
        $_bcs(shape::Tuple, ::Tuple{}) = (shape[1], $_bcs(Base.tail(shape), ())...)
        function $_bcs(shape::Tuple, newshape::Tuple)
            return ($_bcs1(shape[1], newshape[1]), $_bcs(Base.tail(shape), Base.tail(newshape))...)
        end
        $_bcs1(a::Integer, b::Integer) = $_m(a,b)
        $_bcs1(a, b::Integer) = b
        $_bcs1(a::Integer, b) = a
    end
end


import Base: reshape
reshape(A::DoubleInfiniteVector, a::DoubleInfinity, b::Integer...) = _reshape(A, tuple(a,b...))
reshape(A::DoubleInfiniteVector, a::Integer, b::DoubleInfinity, c::Integer...) = _reshape(A, tuple(a,b,c...))
reshape(A::DoubleInfiniteVector, dims::Tuple{DoubleInfinity,Vararg{Integer}}) = _reshape(A, dims)
reshape(A::DoubleInfiniteVector, dims::Tuple{Integer,DoubleInfinity,Vararg{Integer}}) = _reshape(A, dims)
reshape(A::DoubleInfiniteVector, a::Integer, b::Integer, c::Union{Integer,DoubleInfinity}...) = _reshape(A, tuple(a,b,c...))
reshape(A::DoubleInfiniteVector, dims::Tuple{Integer,Integer,Vararg{Union{Integer,DoubleInfinity}}}) = _reshape(A, dims)

function _reshape(A::CompactInfiniteVector, s::Tuple{Vararg{Union{Integer,DoubleInfinity}}})
    i1 = findfirst(isinf, s)
    i2 = findlast(isinf, s)
    if i1!=i2
        throw(DimensionMismatch("Only one dimenson must be infinite"))
    end
    if !all(x->x==1||isinf(x), s)
        throw(DimensionMismatch("Only one dimenson must be infinite in size, all others 1"))
    end
    CompactInfiniteArray(
        reshape(subarray(A),
        ntuple(k->isinf(s[k]) ? sublength(A) : 1, Val(length(s)))),
        ntuple(k->isinf(s[k]) ? offset(A) : 0, Val(length(s))))
end

"""
    δ(n::Int)

Construct Dirac delts, i.e., a bi-infinite vector with δ(k)=1, if k=n, and δ(k)=0, otherwise.
"""
δ(k::Int) = CompactInfiniteVector([1], k)


CompactInfiniteVector(a::Vector{T}, offset = 0) where {T} = CompactInfiniteVector{T}(a, offset)
CompactInfiniteVector(a::AbstractVector{T}, offset = 0) where {T} = CompactInfiniteVector{T}(collect(a), offset)

@inline _offset(vec::CompactInfiniteArray) = vec.offset
@inline offset(vec::CompactInfiniteVector) = vec.offset[1]
@inline offset(vec::CompactInfiniteArray) = vec.offset


@inline setindex!(vec::CompactInfiniteArray{T,N}, k::Vararg{Integer,N}) where {T,N} =
    all( (k .< _offset(vec)) .| (k .>= _offset(vec).+subsize(vec))) ? error("Not possible") : setindex!(subarray(vec), _mapindex(vec, k))


"""
    shift!(vec::InfiniteArrays, m::Int)

In-place shifting of vector. (not always possible)
"""
shift!(vector::CompactInfiniteArray{T,N}, k::Integer) where {T,N} = shift!(vector, ntuple(i->k,Val(N)))
shift!(vector::CompactInfiniteArray{T,N}, k::NTuple{N,Integer}) where {T,N} = (vector.offset = _offset(vector).+k; vector)

"""
    reverse!(vec::CompactInfiniteVector)

In-place time reversel: `vec(-k)`
"""
reverse!(vec::CompactInfiniteArray) = (reverse!(subarray(vec)); vec.offset = .-(_lastindex(vec));vec)



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
    subarray       ::  SVector{L,T}

    FixedInfiniteVector(a::SVector{L,T}, ::Val{OFS}) where {L,T,OFS} = new{L,OFS,T}(a)
end

# type unstable constructor
"""
    FixedInfiniteVector(a::AbstractVector{T}, offset = 0)

Construct a FixedInfiniteVector with indices starting at offset.
"""
FixedInfiniteVector(a::AbstractVector, offset::Integer = 0) = FixedInfiniteVector(length(a)==0 ? SVector{0,eltype(a)}() : SVector(a...), Val{offset}())
FixedInfiniteVector(a::AbstractVector, offset::Tuple{Integer}) =
    FixedInfiniteVector(a,offset[1])

# type unstable constructor
FixedInfiniteVector(vec::CompactInfiniteVector) = FixedInfiniteVector(subarray(vec), _offset(vec))

_offset(vec::FixedInfiniteVector{L,OFS}) where {L,OFS} = (OFS,)
offset(vec::FixedInfiniteVector{L,OFS}) where {L,OFS} = OFS

sublength(vec::FixedInfiniteVector{L,OFS,T}) where {L,OFS,T} = L

function Broadcast.materialize(bc::Broadcast.Broadcasted{<:InfiniteArrayStyle{<:FixedInfiniteVector}})
    ixs = eachnonzeroindex.(bc.args)
    i1s = _bcfirst(bc.args...)
    ins = _bclast(bc.args...)
    sarrays = map(x->_resize(ins.-i1s.+1,i1s,__first(x),subarray(x)),bc.args)
    FixedInfiniteVector(SVector(Broadcast.broadcast(bc.f,  sarrays...)...), Val(i1s[1]))
end

for COMPACTVECTOR in (:CompactInfiniteVector,:FixedInfiniteVector)
    @eval shift(vector::$COMPACTVECTOR, k::Integer) = $COMPACTVECTOR(copy(subarray(vector)), offset(vector)+k)
    @eval reverse(vec::$COMPACTVECTOR) = $COMPACTVECTOR(reverse(subarray(vec), 1), 0 .-(_lastindex(vec)))
    @eval downsample(vec::$COMPACTVECTOR, m::Integer) =
        downsample(vec, ntuple(k->m,Val(ndims(vec))))
    @eval upsample(vec::$COMPACTVECTOR, m::Integer) =
        upsample(vec, ntuple(k->m,Val(ndims(vec))))
    @eval function alternating_flip(vec::$COMPACTVECTOR, pivot = 1)
        hflip = similar(subarray(vec))
        # The first element of hflip has index pivot-_lastindex(f).
        # Whether or not we need to flip its sign depend on the parity of this index:
        isodd(pivot-_lastindex(vec)[1]) ? t = -1 : t = 1
        hflip[1:2:end] =  t * subarray(vec)[end:-2:1]
        hflip[2:2:end] = -t * subarray(vec)[end-1:-2:1]
        $COMPACTVECTOR(hflip, pivot .- _lastindex(vec))
    end
    @eval function alternating(vec::$COMPACTVECTOR)
        halt = similar(subarray(vec))
        t = (-1)^_firstindex(vec)[1]
        for k in eachindex(subarray(vec))
            halt[k] = t * subarray(vec)[k]
            t = -t
        end
        $COMPACTVECTOR(halt, _firstindex(vec))
    end
    @eval support(vec::$COMPACTVECTOR, j::Int, k::Int) = (1/(1<<j)*(support(vec)[1][1]+k), 1/(1<<j)*(support(vec)[1][2]+k))
    @eval function inv(a::$COMPACTVECTOR{T}, q::Int, tol = sqrt(eps(T)); K=sublength(a), maximum=Inf) where {T}
        if q == 1
            return inv(a)
        end
        n = sublength(a)
        Q = iseven(n) ? n>>1 : (n-1) >>1

        # Symmetrise such that I ≈ -n/2..n/2
        sym_shift = -Q-offset(a)
        b = shift(a, sym_shift)
        I = -Q:Q
        L = fld(K+Q,q)

        Ii = -L:L
        Ij = -K:K
        A = zeros(T, length(Ii), length(Ij))
        for i in Ii
            for j in Ij
                A[i-Ii[1]+1, j-Ij[1]+1] = b[q*i-j]
            end
        end
        M = norm(a, Inf)
        A .= A./ M
        # determine rhs
        e = δ(0)[Ii]
        result = pinv(A,tol)*e
        if norm(result, Inf) >= maximum
            error("Can not find compact dual (maximum is $(norm(result, Inf)) (not below $maximum))")
        end
        res = norm(A*result-e)
        if res>tol
            error("Can not find compact dual (residual is $(res), not $tol)")
        end
        # solve and shift back
        $COMPACTVECTOR(result ./ M, sym_shift+first(Ij))
    end
end

for COMPACTARRAY in (:CompactInfiniteArray,:FixedInfiniteVector)
    for fun in (:copy, :similar)
        @eval Base.$fun(vec::$COMPACTARRAY) = $COMPACTARRAY(Base.$fun(subarray(vec)), _offset(vec))
    end
    @eval eachnonzeroindex(vec::$COMPACTARRAY) = eachindex(IndexCartesian(),subarray(vec)) .+ CartesianIndex(_offset(vec).-1)
    @eval shift(vector::$COMPACTARRAY, k::Tuple{Vararg{Integer}}) = $COMPACTARRAY(copy(subarray(vector)), offset(vector).+k)
    @eval LinearAlgebra.norm(vec::$COMPACTARRAY, p::Real=2) = norm(subarray(vec), p)
    @eval Base.minimum(vec::$COMPACTARRAY) = minimum(subarray(vec))
    @eval downsample(vec::$COMPACTARRAY, m::Tuple{Vararg{Integer}}) =
        $COMPACTARRAY(subarray(vec)[ntuple(k->((mod(m[k]-_offset(vec)[k], m[k])+1):m[k]:size(subarray(vec),k)),Val(ndims(vec)))...], cld.(_offset(vec), m))
    @eval function upsample(vec::$COMPACTARRAY, m::Tuple{Vararg{Integer}})
        v = similar(vec, m.*(subsize(vec).-1).+1)
        fill!(v, 0)
        v[ntuple(k->(1:m[k]:size(v,k)),Val(ndims(vec)))...] .= subarray(vec)
        $COMPACTARRAY(v, m.*_offset(vec))
    end
    @eval @inline _mapindex(vec::$COMPACTARRAY, k) = (k) .- (_offset(vec)) .+ 1
    @eval @inline i_mapindex(vec::$COMPACTARRAY, l) = (l) .+ (_offset(vec)) .- 1
    @eval @inline _firstindex(vec::$COMPACTARRAY) = i_mapindex(vec, ntuple(k->1,Val(ndims(vec))))
    @eval @inline _lastindex(vec::$COMPACTARRAY) = i_mapindex(vec, subsize(vec))
    @eval @inline getindex(vec::$COMPACTARRAY, k::Vararg{Integer}) =
        all(ki<oi || ki>=oi+si for (ki,oi,si) in zip(k,offset(vec),subsize(vec))) ? zero(eltype(vec)) : getindex(subarray(vec), _mapindex(vec, k)...)
    @eval function getindex(vec::$COMPACTARRAY, c::Vararg{UnitRange})
        e = zeros(eltype(vec), length.(c))
        c1 = [first(ci) for ci in c]
        i1 = max.(_offset(vec), c1)
        i2 = min.(_offset(vec) .+ subsize(vec).-1, [last(ci) for ci in c])
        e1 = max.(1, 1 .-c1.+_offset(vec))
        e[ntuple(k->((0:(i2[k]-i1[k])).+e1[k]),Val(ndims(vec)))...] .= subarray(vec)[ntuple(k->((i1[k]:i2[k]).+(-_offset(vec)[k]+1)),Val(ndims(vec)))...]
        e
    end
    @eval function evenpart(vec::$COMPACTARRAY)
        I = Iterators.product(map((f,l)->nexteven(f):2:l,_firstindex(vec),_lastindex(vec)))
        $COMPACTARRAY([vec[j...] for j in I][1], div.(nexteven.(_firstindex(vec)),2))
    end
    @eval function oddpart(vec::$COMPACTARRAY)
        I = Iterators.product(map((f,l)->nextodd(f):2:l,_firstindex(vec),_lastindex(vec)))
            $COMPACTARRAY([vec[j...] for j in I][1], div.(previouseven.(_firstindex(vec)),2))
    end
    @eval support(vec::$COMPACTARRAY) = ntuple(k->_offset(vec)[k]:_offset(vec)[k]+subsize(vec)[k]-1,Val(ndims(vec)))
    for op in (:*,:+,:-,:/)
        @eval $op(vec::$COMPACTARRAY, a::Number) = $COMPACTARRAY($op(subarray(vec),a), _offset(vec))
        @eval $op(a::Number, vec::$COMPACTARRAY) = $COMPACTARRAY($op(a,subarray(vec)), _offset(vec))
    end
    for op in (:+,:-)
        @eval function $op(vec1::$COMPACTARRAY, vec2::$COMPACTARRAY)
            f1 = min.(_firstindex(vec1), _firstindex(vec2))
            f2 = max.(_lastindex(vec1), _lastindex(vec2))
            a = zeros(promote_type(eltype(vec1),eltype(vec2)),f2.-f1.+1)
            for i in eachnonzeroindex(vec1)
                a[ntuple(k->i[k]-f1[k]+1,Val(ndims(vec1)))...] = vec1[i]
            end
            for i in eachnonzeroindex(vec2)
                ix = ntuple(k->i[k]-f1[k]+1,Val(ndims(vec1)))
                a[ix...] = $op(a[ix...],vec2[i])
            end
            $COMPACTARRAY(a, f1)
        end
    end
end
@eval conv(vec1::Union{CompactInfiniteVector,FixedInfiniteVector}, vec2::Union{CompactInfiniteVector,FixedInfiniteVector}) =
    CompactInfiniteVector(conv(convert(Vector, copy(subarray(vec1))), convert(Vector, copy(subarray(vec2)))), _offset(vec1) .+ _offset(vec2))

@eval conv(vec1::FixedInfiniteVector, vec2::FixedInfiniteVector) =
    FixedInfiniteVector(conv(convert(Vector,subarray(vec1)), convert(Vector,subarray(vec2))), _offset(vec1) + _offset(vec2))



"The first even number greater than or equal to n."
nexteven(n) = isodd(n) ? n+1 : n

"The last even number, smaller than or equal to n."
previouseven(n) = isodd(n) ? n-1 : n

"The first odd number greater than or equal to n."
nextodd(n) = isodd(n) ? n : n+1

export CompactPeriodicInfiniteVector
"""
    mutable struct CompactPeriodicInfiniteVector{T} <: AbstractFiniteNZInfiniteVector{T}

A CompactPeriodicInfiniteVector contains a sequence of nonzero elements starting at a given offset and has a given period.

    CompactPeriodicInfiniteVector(a::AbstractVector{T}, period, offset = 0)

Construct a CompactPeriodicInfiniteVector with indices `a` of a period `period` starting at `offset`.

# Examples
```jldoctest; setup = :(using InfiniteVectors)
julia> CompactPeriodicInfiniteVector(1:3,,5,-1)
CompactInfiniteVector{Int64} with indices ℤ:
[  …, 2, 3, 0, 0, 1, 2, 3, 0, 0, 1, …  ]
```
"""
mutable struct CompactPeriodicInfiniteArray{T,N,A<:Array{T,N}} <: AbstractFiniteNZInfiniteArray{T,N}
    subarray  ::  A
    period  ::  NTuple{N,Int}
    offset  ::  NTuple{N,Int}

    function CompactPeriodicInfiniteArray{T,1}(a, period::Int, offset::Int) where {T}
        b = zeros(T, min(period,length(a)))
        i = 1
        for ai in a
            b[i] += ai
            i+=1
            if i>period
                i=1
            end
        end
        new{T,1,typeof(b)}(b, tuple(period), tuple(offset))
    end

    function CompactPeriodicInfiniteArray{T,N}(a, period::NTuple{N,Int}, offset::NTuple{N,Int}) where {T,N}
        b = zeros(T, min.(period,size(a)))
        for (i,ai) in zip(CartesianIndices(axes(a)),a)
            b[CartesianIndex( mod.(i.I .- 1,period).+1)] = ai
        end
        new{T,N,typeof(b)}(b, period, offset)
    end
end

CompactPeriodicInfiniteArray(a,period,offset) =
    CompactPeriodicInfiniteArray{eltype(a),ndims(a)}(a, period, offset)
const CompactPeriodicInfiniteVector{T} = CompactPeriodicInfiniteArray{T,1}
CompactPeriodicInfiniteVector(a,period,offset) =
    CompactPeriodicInfiniteArray{eltype(a),1}(a, period, offset)
# type unstable constructor
"""
    CompactPeriodicInfiniteVector(a::AbstractVector{T}, period, offset = 0)

Construct a CompactPeriodicInfiniteVector with indices with period starting at offset.
"""
CompactPeriodicInfiniteVector(a::Vector{T}, period::Int, offset::Int = 0) where {T} = CompactPeriodicInfiniteVector{T}(a, period, offset)
CompactPeriodicInfiniteVector(a::AbstractVector{T}, period::Int, offset::Int = 0) where {T} = CompactPeriodicInfiniteVector{T}(collect(a), period, offset)
CompactPeriodicInfiniteVector(a::CompactInfiniteVector{T}, period::Int) where {T} = CompactPeriodicInfiniteVector{T}(copy(subarray(a)), period, offset(a))

copy(vec::CompactPeriodicInfiniteArray) = CompactPeriodicInfiniteArray(copy(subarray(vec)), _period(vec), _offset(vec))

function Broadcast.materialize(bc::Broadcast.Broadcasted{<:InfiniteArrayStyle{<:CompactPeriodicInfiniteVector}})
    p1 = period(bc.args[1])
    if all(p1==period(a) for a in Base.tail(bc.args))
        ixs = eachnonzeroindex.(bc.args)
        i1s = _bcfirst(bc.args...)
        ins = _bclast(bc.args...)
        sarrays = map(x->_resize(ins.-i1s.+1,i1s,__first(x),subarray(x)),bc.args)
        CompactPeriodicInfiniteVector(Broadcast.broadcast(bc.f,  sarrays...),tuple(p1),i1s)
    else
        p = lcm(period.(bc.args)...)
        arrays = [SubArray(a,(0:p-1,)) for a in bc.args]
        PeriodicInfiniteVector(Broadcast.broadcast(bc.f,  arrays...))
    end
end

@inline _offset(vec::CompactPeriodicInfiniteArray) = vec.offset
@inline offset(vec::CompactPeriodicInfiniteVector) = vec.offset[1]
@inline offset(vec::CompactPeriodicInfiniteArray) = vec.offset
@inline period(vec::CompactPeriodicInfiniteVector) = vec.period[1]
@inline _period(vec::CompactPeriodicInfiniteArray) = vec.period
@inline period(vec::CompactPeriodicInfiniteArray) = vec.period

@inline function setindex!(vec::CompactPeriodicInfiniteArray, k::Vararg{Integer})
    error("Not implemented")
end
shift!(vector::CompactPeriodicInfiniteArray{T,N}, k::Integer) where {T,N} =
    shift!(vector,ntuple(i->k,Val(N)))
shift!(vector::CompactPeriodicInfiniteArray{T,N}, k::NTuple{N,Integer}) where {N,T} =
    (vector.offset = _offset(vector).+k; vector)
shift(vec::CompactPeriodicInfiniteArray, k) = shift!(copy(vec), k)


reverse!(vec::CompactPeriodicInfiniteArray) = (reverse!(subarray(vec)); vec.offset = .-(_lastindex(vec));vec)
reverse(vec::CompactPeriodicInfiniteArray) = reverse!(copy(vec))
function alternating(vec::CompactPeriodicInfiniteVector)
    if iseven(period(vec))
        halt = similar(subarray(vec))
        t = (-1)^_firstindex(vec)[1]
        for k in eachindex(subarray(vec))
            halt[k] = t * subarray(vec)[k]
            t = -t
        end
        CompactPeriodicInfiniteVector(halt, period(vec), offset(vec))
    else
        halt = zeros(eltype(vec),sublength(vec)+period(vec))
        t = (-1)^_firstindex(vec)[1]
        for k in eachindex(subarray(vec))
            halt[k] = t * subarray(vec)[k]
            t = -t
        end
        t = (-1)^_firstindex(vec)[1]
        t = -t
        p = period(vec)
        for k in eachindex(subarray(vec))
            halt[k+p] = t * subarray(vec)[k]
            t = -t
        end
        CompactPeriodicInfiniteVector(halt, 2period(vec), offset(vec))
    end
end

function alternating_flip(vec::CompactPeriodicInfiniteVector, pivot=1)
    if iseven(period(vec))
        hflip = similar(subarray(vec))
        isodd(pivot-_lastindex(vec)[1]) ? t = -1 : t = 1
        hflip[1:2:end] .=  t * subarray(vec)[end:-2:1]
        hflip[2:2:end] .= -t * subarray(vec)[end-1:-2:1]
        CompactPeriodicInfiniteVector(hflip, period(vec), pivot-_lastindex(vec)[1])
    else
        hflip = zeros(eltype(vec), sublength(vec)+period(vec))
        isodd(pivot-_lastindex(vec)[1]) ? t = -1 : t = 1
        hflip[1:2:sublength(vec)] .=  t * subarray(vec)[end:-2:1]
        hflip[2:2:sublength(vec)] .= -t * subarray(vec)[end-1:-2:1]

        hflip[period(vec) .+ (1:2:sublength(vec))] .= -t * subarray(vec)[end:-2:1]
        hflip[period(vec) .+ (2:2:sublength(vec))] .= t * subarray(vec)[end-1:-2:1]

        CompactPeriodicInfiniteVector(hflip, 2period(vec), pivot-_lastindex(vec)[1])
    end
end

@inline function getindex(vec::CompactPeriodicInfiniteArray{T,N}, k::Vararg{Integer,N}) where {T,N}
    fldi = fld.(k.-vec.offset,vec.period)
    k = k .- fldi .* vec.period
    all((k .< _offset(vec)) .| (k .>= _offset(vec).+subsize(vec))) ? zero(eltype(vec)) : getindex(subarray(vec), _mapindex(vec, k)...)
end

@inline _mapindex(vec::CompactPeriodicInfiniteArray{T,1}, k::Integer) where {T,N} =
    _mapindex(vec,tuple(k))
@inline i_mapindex(vec::CompactPeriodicInfiniteArray{T,1}, k::Integer) where {T,N} =
    i_mapindex(vec,tuple(k))
@inline _mapindex(vec::CompactPeriodicInfiniteArray{T,N}, k::NTuple{N,Integer}) where {T,N} =
    k .- _offset(vec) .+ 1
@inline i_mapindex(vec::CompactPeriodicInfiniteArray{T,N}, k::NTuple{N,Integer}) where {T,N} =
    k .+ _offset(vec) .- 1
@inline _firstindex(vec::CompactPeriodicInfiniteArray{T,N}) where {T,N} = i_mapindex(vec, ntuple(k->1,Val(N)))
@inline _lastindex(vec::CompactPeriodicInfiniteArray{T,N}) where {T,N} = i_mapindex(vec, subsize(vec))
conv(vec1::CompactPeriodicInfiniteVector, vec2::CompactPeriodicInfiniteVector) =
    circconv(vec1,vec2)
circconv(vec1::CompactPeriodicInfiniteVector, vec2::CompactPeriodicInfiniteVector) =
        CompactPeriodicInfiniteVector(conv(copy(subarray(vec1)), copy(subarray(vec2))),lcm(vec1.period,vec2.period),offset(vec1) + offset(vec2))
eachnonzeroindex(vec::CompactPeriodicInfiniteArray) = eachindex(subarray(vec)) .+ (_offset(vec) .- 1)

upsample(vec::CompactPeriodicInfiniteArray{T,N}, i::Integer) where {T,N} =
    upsample(vec, ntuple(k->i,Val(N)))
function upsample(vec::CompactPeriodicInfiniteArray{T,N}, m::NTuple{N,Integer}) where {T,N}
    v = similar(vec, m.*(subsize(vec).-1).+1)
    fill!(v, 0)
    v[ntuple(k->1:m[k]:size(v,k),Val(N))...] .= subarray(vec)
    CompactPeriodicInfiniteArray(v, m.*_period(vec), m.*_offset(vec))
end

downsample(vec::CompactPeriodicInfiniteArray{T,N}, i::Integer) where {T,N} =
    downsample(vec, ntuple(k->i,Val(N)))
function downsample(vec::CompactPeriodicInfiniteArray{T,N}, m::NTuple{N,Integer}) where {T,N}
    if  all(x->x==0,mod.(_period(vec), m))
        CompactPeriodicInfiniteVector(subarray(vec)[ntuple(k->mod(m[k]-_offset(vec)[k], m[k])+1:m[k]:subsize(vec)[k],Val(N))...],
            div.(_period(vec),m), cld.(_offset(vec), m))
    else
        downsample(PeriodicInfiniteVector(vec),m)
    end
end

function evenpart(vec::CompactPeriodicInfiniteArray)
    if all(iseven.(_period(vec)))
        I = Iterators.product(map((x,y) -> nexteven(x):2:y ,_firstindex(vec) ,_lastindex(vec)))
        CompactPeriodicInfiniteArray([vec[i...] for i in I ][1], _period(vec).>>1, div.(nexteven.(_firstindex(vec)),2))
    else
        I = Iterators.product(map((x,y,z) -> nexteven(x):2:(z+y) ,_firstindex(vec) ,_lastindex(vec),_period(vec)))
        CompactPeriodicInfiniteArray([vec[i...] for i in I ][1], _period(vec), div.(nexteven.(_firstindex(vec)),2))
        # CompactPeriodicInfiniteVector([vec[j] for j in nexteven(_firstindex(vec)):2:period(vec)+_lastindex(vec)], period(vec), div(nexteven(_firstindex(vec)),2))
    end
end

function oddpart(vec::CompactPeriodicInfiniteArray)
    if  all(iseven.(_period(vec)))
        I = Iterators.product(map((x,y) -> nextodd(x):2:y ,_firstindex(vec) ,_lastindex(vec)))
        CompactPeriodicInfiniteArray([vec[i...] for i in I ][1], _period(vec).>>1, div.(previouseven.(_firstindex(vec)),2))
    else
        I = Iterators.product(map((x,y,z) -> nextodd(x):2:(y+z) ,_firstindex(vec) ,_lastindex(vec),_period(vec)))
        CompactPeriodicInfiniteArray([vec[i...] for i in I ][1], _period(vec), div.(previouseven.(_firstindex(vec)),2))
    end
end

@doc """
    reverse(vec::DoubleInfiniteVector)

Time-reversel: `vec(-k)`
"""
reverse(vec)

@doc """
    alternating_flip(vec::DoubleInfiniteVector, pivot = 1)

From a given filter \$h(i)\$, compute a new filter satisfying the alternating flip relation,
centered around the given pivot:

``g(k) = (-1)^k h(pivot-k)``
"""
alternating_flip(vec, pivot = 1)

@doc """
    alternating(vec::DoubleInfiniteVector)

Compute a new 'alternating' filter satisfying \$g(k) = (-1)^k h(k)`.\$
"""
alternating(vec)

@doc """
evenpart(vec::DoubleInfiniteVector)

    Return the even part of a sequence `s`, defined by `s_e[k] = s[2k]`.
"""
evenpart(vec)

@doc """
    oddpart(vec::DoubleInfiniteVector)

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
