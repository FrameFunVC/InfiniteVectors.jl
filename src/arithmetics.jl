
conv(v1::PeriodicInfiniteVector, v2::CompactInfiniteVector) =
    conv(v2, v1)

function conv(v1::CompactInfiniteVector, v2::PeriodicInfiniteVector)
    r = conv(subvector(v1), subvector(v2))
    p = period(v2)
    v = zeros(eltype(r), p)
    for i in length(v)+1:length(r)
        r[mod(i-1, p)+1] += r[i]
    end
    PeriodicInfiniteVector(circshift(view(r, 1:length(v)), offset(v1)))
end
