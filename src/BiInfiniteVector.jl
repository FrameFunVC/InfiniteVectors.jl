
"""
Any subtype of InfiniteVector is a bi-infinite list of values, with one value for each integer in Z.

They are not iterable
"""
abstract type BiInfiniteVector{T} <: InfiniteVector{T} end
# AbstractVector interface
Base.eachindex(::BiInfiniteVector) = Integers()
Base.axes1(::BiInfiniteVector) = Integers()
Base.axes(::BiInfiniteVector) = (Integers(),)
Base.checkbounds(::Type{Bool}, A::BiInfiniteVector, I...) = all(in.(Ref(eachindex(A)), I))
Base.BroadcastStyle(::Type{T}) where {T<:InfiniteVector}= Broadcast.ArrayStyle{T}()
Base.print_array(io::IO, X::BiInfiniteVector) = Base.show_vector(io, X)
Base.LinearIndices(::NTuple{N,Integers}) where N = Integers()

function Base.show_vector(io::IO, v::BiInfiniteVector, opn='[', cls=']')
    print(io, Base.typeinfo_prefix(io, v))
    # directly or indirectly, the context now knows about eltype(v)
    io = IOContext(io, :typeinfo => eltype(v), :compact => get(io, :compact, true))
    Base.show_delim_array(io, v, opn* "  …, ", ",", ", …  "*cls, false, 0, 9)
end

"""
    transpose(vec::BiInfiniteVector)

Returns the time-reversed vector `vec(-k)`
"""
transpose(vec::BiInfiniteVector) = reverse(vec)

"""
    adjoint(vec::BiInfiniteVector)

Returns the paraconjugate \$\\overline{vec(-k)}\$`
"""
adjoint(vec::BiInfiniteVector) = reverse(conj.(vec))

getindex(vec::BiInfiniteVector, range) = eltype(vec)[vec[k] for k in range]



"""
    fouriertransform(s::BiInfiniteVector, ω)

The Fourier transform of a sequence is defined by \$S(ω) = \\sum_{k\\in ℤ} s_k e^{-i ω k}\$. It is like
the Z transform with \$z = e^{i ω}\$. The Fourier transform is a `2π`-periodic continuous function
of `ω`.
"""
fouriertransform(s::InfiniteVector, ω) = ztransform(s, exp(im*ω))

"""
    ZTransform{S<:BiInfiniteVector} <: Function
ZTransform is a wrapper type for the `ztransform` function.
"""
struct ZTransform{S<:BiInfiniteVector} <: Function
    vec     ::  S
end

(f::ZTransform)(z) = ztransform(z.vec, z)

struct Convolution{T,S1<:BiInfiniteVector,S2<:BiInfiniteVector} <: BiInfiniteVector{T}
    vec1  ::  S1
    vec2  ::  S2
end

Base.getindex(::Convolution) = throw(MethodError("Not implemented"))

conv(vec1::BiInfiniteVector, vec2::BiInfiniteVector) =
    Convolution{promote_type(eltype(vec1),eltype(vec2)), typeof(vec1), typeof(vec2)}(vec1, vec2)

ztransform(vec::Convolution, z) = ztransform(vec.vec1, z) * ztransform(vec.vec2, z)

"""
    leastsquares_inv(v::BiInfiniteVector, m::Int)

The inverse filter
``f = v*[([v*v]_{↓m})^{-1}]_{↑m}``

returns by applying on b
``v*[[f*b]_{↓m}]_{↑m}``
the least squares approximation of `b` in the range of `v`.
"""
leastsquares_inv(v::BiInfiniteVector, m::Int) = v*upsample(inv(downsample(v*v, m)), m)


# function getindex(s::Convolution, k)
#     #TODO: implement
# end
