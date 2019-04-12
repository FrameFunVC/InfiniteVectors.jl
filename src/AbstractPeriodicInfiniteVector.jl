abstract type AbstractPeriodicInfiniteVector{T} <: AbstractDoubleInfiniteVector{T} end


@inline period(vec::AbstractPeriodicInfiniteVector) = length(subvector(vec))

"The vector that goes from 0 to P-1, where P is the period of `vec`"
@inline subvector(vec::AbstractPeriodicInfiniteVector) = vec.subvec

mutable struct PeriodicInfiniteVector{T} <: AbstractPeriodicInfiniteVector{T}
    subvec      :: Vector{T}
end

@inline getindex(vec::PeriodicInfiniteVector, i) = getindex(subvector(vec), mod.(i, period(vec)))
@inline setindex!(vec::PeriodicInfiniteVector, i) = setindex!(subvector(vec), mod.(i, period(vec)))

conv(vec1::PeriodicInfiniteVector, vec2::PeriodicInfiniteVector) =
    period(vec1)==period(vec2) ?
        PeriodicInfiniteVector(periodic_conv(subvector(vec1),subvector(vec2))) :
        Convolution(vec2,vec2)

shift(vector::PeriodicInfiniteVector, k::Int) = PeriodicInfiniteVector(circshift(subvector(vec), k))
shift!(vector::PeriodicInfiniteVector, k::Int) = circshift!(subvector(vec), k)

function periodic_conv(u::StridedVector{T}, v::StridedVector{T}) where T<:BLAS.BlasFloat
    nu = length(u)
    nv = length(v)
    if nu != mv
        error("Periodic convolution not posssible")
    end
    n = nu
    np2 = n > 1024 ? nextprod([2,3,5], n) : nextpow(2, n)
    upad = [u; zeros(T, np2 - nu)]
    vpad = [v; zeros(T, np2 - nv)]
    if T <: Real
        p = plan_rfft(upad)
        y = irfft((p*upad).*(p*vpad), np2)
    else
        p = plan_fft!(upad)
        y = ifft!((p*upad).*(p*vpad))
    end
    return y[1:n]
end
