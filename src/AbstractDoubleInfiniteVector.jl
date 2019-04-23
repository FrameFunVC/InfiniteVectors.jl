
"""
Any subtype of InfiniteVector is a bi-infinite list of values, with one value for each integer in Z.

They are not iterable
"""
abstract type AbstractDoubleInfiniteVector{T} <: InfiniteVector{T} end
# AbstractVector interface
Base.eachindex(::AbstractDoubleInfiniteVector) = Integers()
Base.axes1(::AbstractDoubleInfiniteVector) = Integers()
Base.axes(::AbstractDoubleInfiniteVector) = (Integers(),)
Base.checkbounds(::Type{Bool}, A::AbstractDoubleInfiniteVector, I...) = all(in.(Ref(eachindex(A)), I))
Base.BroadcastStyle(::Type{T}) where {T<:InfiniteVector}= Broadcast.ArrayStyle{T}()
Base.print_array(io::IO, X::AbstractDoubleInfiniteVector) = Base.show_vector(io, X)

function Base.show_vector(io::IO, v::AbstractDoubleInfiniteVector, opn='[', cls=']')
    print(io, Base.typeinfo_prefix(io, v))
    # directly or indirectly, the context now knows about eltype(v)
    io = IOContext(io, :typeinfo => eltype(v), :compact => get(io, :compact, true))
    Base.show_delim_array(io, v, opn* "  …, ", ",", ", …  "*cls, false, 0, 9)
end

# transpose is time reverse
transpose(vec::AbstractDoubleInfiniteVector) = reverse(vec)
# adjoint is
adjoint(vec::AbstractDoubleInfiniteVector) = reverse(conj.(vec))

getindex(vec::AbstractDoubleInfiniteVector, range) = eltype(vec)[vec[k] for k in range]

function upsample(vec::AbstractDoubleInfiniteVector, m::Int) end

function downsample(vec::AbstractDoubleInfiniteVector, m::Int) end

"The j-th discrete moment of a sequence is defined as `\\sum_k h_k k^j`."
function moment(vec::AbstractDoubleInfiniteVector, j) end

"""
The Z transform of a sequence is a continuous function of `z`, defined by
`S(z) = \\sum_k s_k z^{-k}`.
"""
function ztransform(vec::AbstractDoubleInfiniteVector, z) end

"""
The Fourier transform of a sequence is defined by `S(ω) = \\sum_k s_k e^{-i ω k}`. It is like
the Z transform with `z = e^{i ω}`. The Fourier transform is a `2π`-periodic continuous function
of `ω`.
"""
fouriertransform(s::AbstractDoubleInfiniteVector, ω) = ztransform(s, exp(im*ω))

"ZTransform is a wrapper type for the `ztransform` function."
struct ZTransform{S<:AbstractDoubleInfiniteVector} <: Function
    vec     ::  S
end

(f::ZTransform)(z) = ztransform(z.vec, z)

struct Convolution{T,S1<:AbstractDoubleInfiniteVector,S2<:AbstractDoubleInfiniteVector} <: AbstractDoubleInfiniteVector{T}
    vec1  ::  S1
    vec2  ::  S2
end

conv(vec1::AbstractDoubleInfiniteVector, vec2::AbstractDoubleInfiniteVector) =
    Convolution{promote_type(eltype(vec1),eltype(vec2)), typeof(vec1), typeof(vec2)}(vec1, vec2)

ztransform(vec::Convolution, z) = ztransform(vec.vec1, z) * ztransform(vec.vec2, z)

"""
The inverse filter `f` such that for some `b`

    v*[[f*b]↓m]↑m

is a least squares approximation of `b` in the range of `v`.
It is given by

    f = v*[([v*v]↓m)^{-1}]↑m
"""
leastsquares_inv(v::AbstractDoubleInfiniteVector, m::Int) = v*upsample(inv(downsample(v*v, m)), m)


# function getindex(s::Convolution, k)
#     #TODO: implement
# end
