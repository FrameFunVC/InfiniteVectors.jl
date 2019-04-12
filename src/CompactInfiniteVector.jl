# compactsequence.jl


#########################################
# CompactInfiniteVector
#########################################

"""
A CompactInfiniteVector is compactly supported sequence of a certain length starting at a
given offset index.
"""
mutable struct CompactInfiniteVector{T} <: AbstractFiniteNZInfiniteVector{T}
    subvec  ::  Vector{T}
    offset  ::  Int

    CompactInfiniteVector{T}(a, offset) where T = new(a, offset)
end

δ(k::Int) = CompactInfiniteVector([1], k)
CompactInfiniteVector(a::Vector{T}, offset = 0) where {T} = CompactInfiniteVector{T}(a, offset)
CompactInfiniteVector(a::AbstractVector{T}, offset = 0) where {T} = CompactInfiniteVector{T}(collect(a), offset)

function copy(bc::Base.Broadcast.Broadcasted{<:Base.Broadcast.ArrayStyle{<:CompactInfiniteVector},<:Tuple,F,Tuple{T}} where {F,T<:CompactInfiniteVector})
    vc = bc.args[1]
    CompactInfiniteVector((bc.f).(subvector(vc)), offset(vc))
end

copy(vec::CompactInfiniteVector) = CompactInfiniteVector(copy(subvector(vec)), offset(vec))

for op in (:(==), :≈,)
    @eval $op(vec1::CompactInfiniteVector, vec2::CompactInfiniteVector) =
        $op(subvector(vec1), subvector(vec2)) && (offset(vec1)==offset(vec2))
end

@inline subvector(vec::CompactInfiniteVector) = vec.subvec
@inline sublength(vec::InfiniteVector) = length(subvector(vec))
@inline offset(vec::CompactInfiniteVector) = vec.offset
eachnonzeroindex(vec::CompactInfiniteVector) = (1:sublength(vec)) .+ offset(vec)


shift(vector::CompactInfiniteVector, k::Int) = CompactInfiniteVector(copy(subvector(vector)), offset(vector)+k)
shift!(vector::CompactInfiniteVector, k::Int) = (vector.offset = offset(vector)+k; vector)

reverse(vec::CompactInfiniteVector) = CompactInfiniteVector(reverse(subvector(vec), 1), -_lastindex(vec))
reverse!(vec::CompactInfiniteVector) = (reverse!(subvector(vec)); vec.offset = -_lastindex(vec);vec)

downsample(vec::CompactInfiniteVector, m::Int) =
    CompactInfiniteVector(subvector(vec)[mod(m-offset(vec), m)+1:m:end], cld(offset(vec), m))

function upsample(vec::CompactInfiniteVector, m::Int)
    v = similar(vec, m*(sublength(vec)-1)+1)
    fill!(v, 0)
    v[1:m:end] .= subvector(vec)
    CompactInfiniteVector(v, m*offset(vec))
end


conv(vec1::CompactInfiniteVector, vec2::CompactInfiniteVector) =
    CompactInfiniteVector(conv(subvector(vec1), subvector(vec2)), offset(vec1) + offset(vec2))



@inline _mapindex(vec::CompactInfiniteVector, k) = k - offset(vec) + 1
@inline i_mapindex(vec::CompactInfiniteVector, l) = l + offset(vec) - 1
@inline _firstindex(vec::CompactInfiniteVector) = i_mapindex(vec, 1)
@inline _lastindex(vec::CompactInfiniteVector) = i_mapindex(vec, sublength(vec))

@inline getindex(vec::CompactInfiniteVector, k::Int) =
    (k < offset(vec)) || (k >= offset(vec)+sublength(vec)) ? zero(eltype(vec)) : getindex(subvector(vec), _mapindex(vec, k))
@inline setindex!(vec::CompactInfiniteVector, k::Int) =
    (k < offset(vec)) || (k >= offset(vec)+sublength(vec)) ? error("Not possible") : setindex!(subvector(vec), _mapindex(vec, k))

function getindex(vec::CompactInfiniteVector, c::UnitRange)
    e = zeros(eltype(vec), length(c))
    i1 = max(offset(vec), c[1])
    i2 = min(offset(vec) + sublength(vec)-1, c[end])
    e1 = max(1, 1-c[1]+offset(vec))
    # e2 = min(length(e),i2-i1+1-c[1]+offset(vec))
    e[(0:(i2-i1)) .+ e1] .= subvector(vec)[(i1:i2) .+ (-offset(vec)+1)]
    e
end

"""
From a given filter h_i, compute a new filter satisfying the alternating flip relation,
centered around the given pivot:

`g_k = (-1)^k h_{pivot-k}`

The default pivot is 1.
"""
function alternating_flip(vec::CompactInfiniteVector, pivot = 1)
    hflip = similar(subvector(vec))
    # The first element of hflip has index pivot-_lastindex(f).
    # Whether or not we need to flip its sign depend on the parity of this index:
    isodd(pivot-_lastindex(vec)) ? t = -1 : t = 1
    hflip[1:2:end] =  t * subvector(vec)[end:-2:1]
    hflip[2:2:end] = -t * subvector(vec)[end-1:-2:1]
    CompactInfiniteVector(hflip, pivot - _lastindex(vec))
end

"Compute a new 'alternating' filter satisfying `g_k = (-1)^k h_k`."
function alternating(vec::CompactInfiniteVector)
    halt = similar(subvector(vec))
    t = (-1)^_firstindex(vec)
    for k in eachindex(subvector(vec))
        halt[k] = t * subvector(vec)[k]
        t = -t
    end
    CompactInfiniteVector(halt, _firstindex(vec))
end


"The first even number greater than or equal to n."
nexteven(n) = isodd(n) ? n+1 : n

"The last even number, smaller than or equal to n."
previouseven(n) = isodd(n) ? n-1 : n

"The first odd number greater than or equal to n."
nextodd(n) = isodd(n) ? n : n+1

"Return the even part of a sequence `s`, defined by `s_e[k] = s[2k]`."
evenpart(vec::CompactInfiniteVector) =
    CompactInfiniteVector(eltype(vec)[vec[j] for j in nexteven(_firstindex(vec)):2:_lastindex(vec)], div(nexteven(_firstindex(vec)),2))

"Return the odd part of a sequence `s`, defined by `s_o[k] = s[2k+1]`."
oddpart(vec::CompactInfiniteVector) =
    CompactInfiniteVector(eltype(vec)[vec[j] for j in nextodd(_firstindex(vec)):2:_lastindex(vec)], div(previouseven(_firstindex(vec)),2))

# conj(vec::CompactInfiniteVector) = CompactInfiniteVector(conj(subvector(vec)), _firstindex(vec))

# moment(vec::CompactInfiniteVector,i) = sum(s[k]*k^i for k in _firstindex(vec):_lastindex(vec))

support(vec::CompactInfiniteVector) = (offset(vec), offset(vec) + sublength(vec) - 1)

support(vec::CompactInfiniteVector, j::Int, k::Int) = (1/(1<<j)*(support(vec)[1]+k), 1/(1<<j)*(support(vec)[2]+k))

sum(vec::CompactInfiniteVector) = sum(subvector(vec))

widen(vec::CompactInfiniteVector) = CompactInfiniteVector(widen(subvector(vec)), offset(vec))
widen(a::Array{T,N}) where {T,N} = convert(Array{widen(T),N}, a)

# for op in (:+, :-, :/, :*)
#     @eval ($op)(x::Number, vec::CompactInfiniteVector) = CompactInfiniteVector(($op)(subvector(vec),x),offset(vec))
#     @eval ($op)(vec::CompactInfiniteVector, x::Number) = ($op)(x, s)
# end

#########################################
# FixedInfiniteVector
#########################################


"""
A FixedInfiniteVector has fixed length `L` starting at a certain offset `OFS`.
It contains its values in a fixed size array.
"""
struct FixedInfiniteVector{L,OFS,T} <: AbstractFiniteNZInfiniteVector{T}
    a       ::  SVector{L,T}
end

FixedInfiniteVector(a::SVector{L,T}, ::Val{OFS}) where {L,T,OFS} = FixedInfiniteVector{L,OFS,T}(a)

# type unsubvector(vec)ble constructor
FixedInfiniteVector(a::AbstractVector, offset = 0) = FixedInfiniteVector(SVector(a...), Val{offset}())

# type unsubvector(vec)ble constructor
FixedInfiniteVector(vec::CompactInfiniteVector) = FixedInfiniteVector(subvector(vec), offset(vec))

convert(::Type{CompactInfiniteVector{T}}, vec::FixedInfiniteVector{L,OFS,T}) where {L,OFS,T}= CompactInfiniteVector(Array(subvector(vec)), OFS)

sublength(vec::FixedInfiniteVector{L,OFS,T}) where {L,OFS,T} = L

getindex(vec::FixedInfiniteVector{L,OFS,T}, i) where {L,OFS,T} = subvector(vec)[i-OFS+1]

_firstindex(vec::FixedInfiniteVector{L,OFS,T}) where {L,OFS,T} = OFS

_lastindex(vec::FixedInfiniteVector{L,OFS,T}) where {L,OFS,T} = OFS+L-1

shift(vec::FixedInfiniteVector{L,OFS,T}, k::Int) where {L,OFS,T} = FixedInfiniteVector(subvector(vec), Val{OFS+k}())



# for op in (:ctranspose, :evenpart, :oddpart, :alternating_flip, :reverse, :conj, :alternating)
#     @eval $op{L,OFS,T}(vec::FixedInfiniteVector{L,OFS,T}) = FixedInfiniteVector($op(CompactInfiniteVector{T}(vec)))
# end
