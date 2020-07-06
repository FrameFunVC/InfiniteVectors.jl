
# AbstractVector interface
Base.eachindex(::DoubleInfiniteVector) = Integers()
Base.print_array(io::IO, X::DoubleInfiniteVector) = Base.show_vector(io, X)
Base.Broadcast.axistype(a::Integers, b) = Integers()
Base.Broadcast.axistype(::T, ::T) where T<:Integers = Integers()
function Base.show_vector(io::IO, v::DoubleInfiniteVector, opn='[', cls=']')
    print(io, Base.typeinfo_prefix(io, v))
    # directly or indirectly, the context now knows about eltype(v)
    io = IOContext(io, :typeinfo => eltype(v), :compact => get(io, :compact, true))
    Base.show_delim_array(io, v, opn* "  …, ", ",", ", …  "*cls, false, 0, 9)
end

"""
    transpose(vec::DoubleInfiniteVector)

Returns the time-reversed vector `vec(-k)`
"""
transpose(vec::DoubleInfiniteVector) = reverse(vec)

"""
    adjoint(vec::DoubleInfiniteVector)

Returns the paraconjugate \$\\overline{vec(-k)}\$`
"""
adjoint(vec::DoubleInfiniteVector) = reverse(conj.(vec))

getindex(vec::DoubleInfiniteVector, range::AbstractArray) = eltype(vec)[vec[k] for k in range]

"""
    fouriertransform(s::DoubleInfiniteVector, ω)

The Fourier transform of a sequence is defined by \$S(ω) = \\sum_{k\\in ℤ} s_k e^{-i ω k}\$. It is like
the Z transform with \$z = e^{i ω}\$. The Fourier transform is a `2π`-periodic continuous function
of `ω`.
"""
fouriertransform(s::DoubleInfiniteVector, ω) = ztransform(s, exp(im*ω))

"""
    ZTransform{S<:DoubleInfiniteVector} <: Function
ZTransform is a wrapper type for the `ztransform` function.
"""
struct ZTransform{S<:DoubleInfiniteVector} <: Function
    vec     ::  S
end

(f::ZTransform)(z) = ztransform(z.vec, z)

struct Convolution{T,S1<:DoubleInfiniteVector,S2<:DoubleInfiniteVector} <: DoubleInfiniteVector{T}
    vec1  ::  S1
    vec2  ::  S2
end

Base.getindex(::Convolution) = throw(MethodError("Not implemented"))

conv(vec1::DoubleInfiniteVector, vec2::DoubleInfiniteVector) =
    Convolution{promote_type(eltype(vec1),eltype(vec2)), typeof(vec1), typeof(vec2)}(vec1, vec2)

ztransform(vec::Convolution, z) = ztransform(vec.vec1, z) * ztransform(vec.vec2, z)

"""
    leastsquares_inv(v::DoubleInfiniteVector, m::Int)

The inverse filter
``f = v*[([v*v]_{↓m})^{-1}]_{↑m}``

returns by applying on b
``v*[[f*b]_{↓m}]_{↑m}``
the least squares approximation of `b` in the range of `v`.
"""
leastsquares_inv(v::DoubleInfiniteVector, m::Int) = v*upsample(inv(downsample(v*v, m)), m)


# function getindex(s::Convolution, k)
#     #TODO: implement
# end
