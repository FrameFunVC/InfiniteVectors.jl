⊛(a, b) = circconv(a, b)
⋆(a, b) = conv(a, b)
*(a::InfiniteVector, b::InfiniteVector) = conv(a, b)
*(a::AbstractPeriodicInfiniteVector, b::AbstractPeriodicInfiniteVector) = circconv(a, b)
*(a::AbstractPeriodicInfiniteVector, b::InfiniteVector) = circconv(b, a)
*(a::InfiniteVector, b::AbstractPeriodicInfiniteVector) = circconv(a, b)


conv(v1::PeriodicInfiniteVector, v2::CompactInfiniteVector) =
    conv(v2, v1)

function circconv(v1::CompactInfiniteVector, v2::PeriodicInfiniteVector)
    r = conv(subvector(v1), subvector(v2))
    p = period(v2)
    v = zeros(eltype(r), p)
    for i in length(v)+1:length(r)
        r[mod(i-1, p)+1] += r[i]
    end
    PeriodicInfiniteVector(circshift(view(r, 1:length(v)), offset(v1)))
end

function PeriodicInfiniteVector(vec::CompactInfiniteVector, N::Int)
    a = zeros(eltype(subvector(vec)), N)
    for i in eachnonzeroindex(vec)
        a[mod(i, N) + 1] += vec[i]
    end
    PeriodicInfiniteVector(a)
end
