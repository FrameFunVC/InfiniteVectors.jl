"""
    abstract type AbstractPeriodicInfiniteVector{T} <: BiInfiniteVector{T} end

A doubly infinite vector that has a period.
"""
abstract type AbstractPeriodicInfiniteVector{T} <: BiInfiniteVector{T} end

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

function copy(bc::Base.Broadcast.Broadcasted{<:Base.Broadcast.ArrayStyle{<:PeriodicInfiniteVector},<:Tuple,F,Tuple{T}} where {F,T<:PeriodicInfiniteVector})
    vc = bc.args[1]
    PeriodicInfiniteVector((bc.f).(subvector(vc)))
end

copy(vec::PeriodicInfiniteVector) = PeriodicInfiniteVector(copy(subvector(vec)))

for op in (:(==), :≈)
    @eval $op(vec1::PeriodicInfiniteVector, vec2::PeriodicInfiniteVector) =
        $op(subvector(vec1), subvector(vec2))
end

@inline getindex(vec::PeriodicInfiniteVector, i) = getindex(subvector(vec), mod.(i, period(vec)) .+ 1)
@inline setindex!(vec::PeriodicInfiniteVector, i) = setindex!(subvector(vec), mod.(i, period(vec)) .+ 1)

circconv(vec1::PeriodicInfiniteVector, vec2::PeriodicInfiniteVector) =
        PeriodicInfiniteVector(circconv(subvector(vec1),subvector(vec2)))
"""
    shift(vec::InfiniteArrays, m::Int)

The resulting vector `r` satisfies `r(k) = vec(k+m)`
"""
shift(vec::PeriodicInfiniteVector, k::Int) = PeriodicInfiniteVector(circshift(subvector(vec), k))

# Not realy in-place
# shift!(vec::PeriodicInfiniteVector, k::Int) = (circshift!(subvector(vec), copy(subvector(vec)), k); vec)

reverse(vec::PeriodicInfiniteVector) = PeriodicInfiniteVector(circshift(reverse(subvector(vec)), 1))
reverse!(vec::PeriodicInfiniteVector) = (circshift!(subvector(vec), reverse!(copy(subvector(vec))), 1);vec)

"""
    downsample(vec::InfiniteArrays, m::Int)

The resulting vector `r` satisfies `r(k) = vec(mk)`
"""
downsample(vec::PeriodicInfiniteVector, m::Int) = mod(period(vec), m)==0 ?
    PeriodicInfiniteVector(subvector(vec)[1:m:end]) :
    PeriodicInfiniteVector(vec[0:m:div(m,gcd(period(vec),m))*(period(vec))-1])

"""
    upample(vec::InfiniteVector, m::Int)

The resulting vector `r` satisfies `r(k) = vec(k/m)`, if `k` is multple of `m`, otherwise, `r(k)=0`
"""
function upsample(vec::PeriodicInfiniteVector, m::Int)
    v = zeros(eltype(subvector(vec)), m*period(vec))
    v[1:m:end] .= subvector(vec)
    PeriodicInfiniteVector(v)
end

function inv(vec::PeriodicInfiniteVector, m::Int)
    return inv(vec)
end

inv(vec::PeriodicInfiniteVector{T}) where T =
    T <: Real ?
        PeriodicInfiniteVector(irfft(1 ./ rfft(subvector(vec)), period(vec))) :
        PeriodicInfiniteVector( ifft(1 ./  fft(subvector(vec))))



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
