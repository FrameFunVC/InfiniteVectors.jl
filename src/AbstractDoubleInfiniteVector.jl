"""
Any subtype of InfiniteVector is a bi-infinite list of values, with one value for each integer in Z.

They are not iterable
"""
abstract type AbstractDoubleInfiniteVector{T} <: AbstractVector{T} end
# AbstractVector interface
Base.size(::AbstractDoubleInfiniteVector) = (∞,)
Base.length(::AbstractDoubleInfiniteVector) = ∞
Base.eachindex(::AbstractDoubleInfiniteVector) = GeometricSpace{Int}()
Base.axes1(::AbstractDoubleInfiniteVector) = GeometricSpace{Int}()
Base.axes(::AbstractDoubleInfiniteVector) = (GeometricSpace{Int}(),)
Base.checkbounds(::Type{Bool}, A::AbstractDoubleInfiniteVector, I...) = all(in.(Ref(eachindex(A)), I))
Base.in(::GeometricSpace{Int}, ::Colon) = true


# transpose is time reverse
transpose(vec::AbstractDoubleInfiniteVector) = reverse(vec)
# adjoint is
adjoint(vec::AbstractDoubleInfiniteVector) = reverse(conj.(vec))

getindex(vec::AbstractDoubleInfiniteVector, range::AbstractRange) = eltype(vec)[vec[k] for k in range]

Base.print_array(io::IO, X::AbstractDoubleInfiniteVector) = Base.show_vector(io, X)

function Base.show_vector(io::IO, v, opn='[', cls=']')
    print(io, Base.typeinfo_prefix(io, v))
    # directly or indirectly, the context now knows about eltype(v)
    io = IOContext(io, :typeinfo => eltype(v), :compact => get(io, :compact, true))
    Base.show_delim_array(io, v, opn* "  …, ", ",", ", …  "*cls, false, 0, 9)
end

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

conv(vec1::AbstractDoubleInfiniteVector, vec2::AbstractDoubleInfiniteVector) = Convolution{promote_eltype(eltype(vec1),eltype(vec2))}(vec1, vec2)

(*)(vec1::AbstractDoubleInfiniteVector, vec2::AbstractDoubleInfiniteVector) = conv(vec1, vec2)

ztransform(vec::Convolution, z) = ztransform(vec.vec1, z) * ztransform(vec.vec2, z)

# function getindex(s::Convolution, k)
#     #TODO: implement
# end
