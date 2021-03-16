"""
    abstract type AbstractPeriodicInfiniteVector{T} <: DoubleInfiniteVector{T} end

A doubly infinite vector that has a period.
"""
abstract type AbstractPeriodicInfiniteVector{T} <: DoubleInfiniteVector{T} end

"""
    period(vec::AbstractPeriodicInfiniteVector)

The period of the periodic vector.
"""
@inline period(vec::AbstractPeriodicInfiniteVector) = length(subvector(vec))

"""
    subvector(vec::AbstractPeriodicInfiniteVector)

The vector that goes from 0 to P-1, where P is the period of `vec`"
"""
@inline subvector(vec::AbstractPeriodicInfiniteVector) = vec.subvec

"""
    mutable struct PeriodicInfiniteVector{T} <: AbstractPeriodicInfiniteVector{T}

A doubly infinite vector that has a period.

PeriodicInfiniteVector(a::AbstractVector{T})

Construct a PeriodicInfiniteVector.
The first element of a becomes the element at index zero, the second one at index 1.

# Examples
```jldoctest; setup = :(using InfiniteVectors)
julia> PeriodicInfiniteVector(1:3)
PeriodicInfiniteVector{Int64} with indices ℤ:
[  …, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, …  ]
```
"""
mutable struct PeriodicInfiniteVector{T} <: AbstractPeriodicInfiniteVector{T}
    subvec      :: Vector{T}
end

PeriodicInfiniteVector(a::AbstractVector{T}) where {T} = PeriodicInfiniteVector{T}(collect(a))

function Broadcast.materialize(bc::Broadcast.Broadcasted{<:InfiniteArrayStyle{<:PeriodicInfiniteVector}})
    p = lcm(period.(bc.args)...)
    arrays = [SubArray(a,(0:p-1,)) for a in bc.args]
    PeriodicInfiniteVector(Broadcast.broadcast(bc.f,  arrays...))
end

function copy(bc::Base.Broadcast.Broadcasted{<:Base.Broadcast.ArrayStyle{<:PeriodicInfiniteVector},<:Tuple,F,Tuple{T}} where {F,T<:PeriodicInfiniteVector})
    vc = bc.args[1]
    PeriodicInfiniteVector((bc.f).(subvector(vc)))
end

copy(vec::PeriodicInfiniteVector) = PeriodicInfiniteVector(copy(subvector(vec)))

for op in (:(==), :≈)
    @eval $op(vec1::PeriodicInfiniteVector, vec2::PeriodicInfiniteVector) =
        $op(subvector(vec1), subvector(vec2))
end

@inline getindex(vec::PeriodicInfiniteVector, i::Union{Integer,Tuple{Vararg{Integer}}}) =
    getindex(subvector(vec), mod.(i, period(vec)) .+ 1)
@inline setindex!(vec::PeriodicInfiniteVector, i::Union{Integer,Tuple{Vararg{Integer}}}) =
    setindex!(subvector(vec), mod.(i, period(vec)) .+ 1)

circconv(vec1::PeriodicInfiniteVector, vec2::PeriodicInfiniteVector) =
        PeriodicInfiniteVector(circconv(subvector(vec1),subvector(vec2)))
"""
    shift(vec::InfiniteArrays, m::Int)

The resulting vector `r` satisfies `r(k) = vec(k+m)`
"""
shift(vec::PeriodicInfiniteVector, k::Integer) =
    shift(vec,ntuple(i->k,Val(ndims(vec))))
shift(vec::PeriodicInfiniteVector, k::Tuple{Vararg{Integer}}) =
    PeriodicInfiniteVector(circshift(subvector(vec), k[1]))
reverse(vec::PeriodicInfiniteVector) = PeriodicInfiniteVector(circshift(reverse(subvector(vec)), 1))

"""
    downsample(vec::InfiniteArrays, m::Int)

The resulting vector `r` satisfies `r(k) = vec(mk)`
"""
downsample(vec::PeriodicInfiniteVector, m::Tuple{Integer}) =
    downsample(vec,first(m))
downsample(vec::PeriodicInfiniteVector, m::Integer) = mod(period(vec), m)==0 ?
    PeriodicInfiniteVector(subvector(vec)[1:m:end]) :
    PeriodicInfiniteVector(vec[0:m:div(m,gcd(period(vec),m))*(period(vec))-1])

"""
    upample(vec::DoubleInfiniteVector, m::Int)

The resulting vector `r` satisfies `r(k) = vec(k/m)`, if `k` is multple of `m`, otherwise, `r(k)=0`
"""
function upsample(vec::PeriodicInfiniteVector, m::Int)
    v = zeros(eltype(subvector(vec)), m*period(vec))
    v[1:m:end] .= subvector(vec)
    PeriodicInfiniteVector(v)
end

for op in (:*,:+,:-,:/)
    @eval $op(vec::PeriodicInfiniteVector, a::Number) = PeriodicInfiniteVector($op(subvector(vec),a), offset(vec))
    @eval $op(a::Number, vec::PeriodicInfiniteVector) = PeriodicInfiniteVector($op(a,subvector(vec)), offset(vec))
end

function inv(vec::PeriodicInfiniteVector, m::Int)
    return inv(vec)
end

function inv(vec::PeriodicInfiniteVector{T}) where T
    if T <: Real
        a = rfft(subvector(vec))
        any(1 .+ a .≈ 1) && error("Singularity error, inverse not possible. ")
        PeriodicInfiniteVector(irfft(1 ./ a, period(vec)))
    else
        a = fft(subvector(vec))
        any(1 .+  a .≈ 1) && error("Singularity error, inverse not possible. ")
        PeriodicInfiniteVector( ifft(1 ./ a))
    end
end


function circconv(u::StridedVector{T}, v::StridedVector{T}) where T<:BLAS.BlasFloat
    nu = length(u)
    nv = length(v)
    n = lcm(nu, nv)
    upad = repeat(u, div(n, nu))
    vpad = repeat(v, div(n, nv))
    if T <: Real
        p = plan_rfft(upad)
        y = irfft((p*upad).*(p*vpad), n)
    else
        p = plan_fft(upad)
        y = ifft!((p*upad).*(p*vpad))
    end
    y
end

function alternating_flip(vec::PeriodicInfiniteVector, pivot = 1)
    r = similar(subvector(vec))
    for i in 0:length(r)-1
        r[i+1] = (-1)^i*vec[pivot-i]
    end
    PeriodicInfiniteVector(r)
end

function alternating(vec::PeriodicInfiniteVector)
    r = copy(subvector(vec))
    for i in 2:2:length(r)
        r[i] *= -1
    end
    PeriodicInfiniteVector(r)
end

evenpart(vec::PeriodicInfiniteVector) =
    iseven(period(vec)) ?
        PeriodicInfiniteVector(subvector(vec)[1:2:end]) :
        PeriodicInfiniteVector(vcat(subvector(vec)[1:2:end], subvector(vec)[2:2:end]))

oddpart(vec::PeriodicInfiniteVector) =
    iseven(period(vec)) ?
        PeriodicInfiniteVector(subvector(vec)[2:2:end]) :
        PeriodicInfiniteVector(vcat(subvector(vec)[2:2:end], subvector(vec)[1:2:end]))
