abstract type AbstractPeriodicInfiniteVector{T} <: AbstractDoubleInfiniteVector{T} end


@inline period(vec::AbstractPeriodicInfiniteVector) = length(subvector(vec))

"The vector that goes from 0 to P-1, where P is the period of `vec`"
@inline subvector(vec::AbstractPeriodicInfiniteVector) = vec.subvec

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

conv(vec1::PeriodicInfiniteVector, vec2::PeriodicInfiniteVector) =
    period(vec1)==period(vec2) ?
        PeriodicInfiniteVector(periodic_conv(subvector(vec1),subvector(vec2))) :
        Convolution(vec2,vec2)

shift(vec::PeriodicInfiniteVector, k::Int) = PeriodicInfiniteVector(circshift(subvector(vec), k))
shift!(vec::PeriodicInfiniteVector, k::Int) = (circshift!(subvector(vec), copy(subvector(vec)), k); vec)

reverse(vec::PeriodicInfiniteVector) = PeriodicInfiniteVector(circshift(reverse(subvector(vec)), 1))
reverse!(vec::PeriodicInfiniteVector) = (circshift!(subvector(vec), reverse!(copy(subvector(vec))), 1);vec)

downsample(vec::PeriodicInfiniteVector, m::Int) = mod(period(vec), m)==0 ?
    PeriodicInfiniteVector(subvector(vec)[1:m:end]) :
    PeriodicInfiniteVector(vec[0:m:div(m,gcd(period(vec),m))*(period(vec))-1])

function upsample(vec::PeriodicInfiniteVector, m::Int)
    v = zeros(eltype(subvector(vec)), m*period(vec))
    v[1:m:end] .= subvector(vec)
    PeriodicInfiniteVector(v)
end

function periodic_conv(u::StridedVector{T}, v::StridedVector{T}) where T<:BLAS.BlasFloat
    nu = length(u)
    nv = length(v)
    n = lcm(nu, nv)
    upad = repeat(u, div(n, nv))
    vpad = repeat(v, div(n, nu))
    if T <: Real
        p = plan_rfft(upad)
        y = irfft((p*upad).*(p*vpad), n)
    else
        p = plan_fft!(upad)
        y = ifft!((p*upad).*(p*vpad))
    end
    y
end
