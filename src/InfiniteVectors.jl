
module InfiniteVectors

using InfiniteArrays, DSP, LinearAlgebra, FFTW, PGFPlotsX
using StaticArrays

import Base: size, length, getindex, setindex!, similar, adjoint, transpose, sum, *, eltype, :,
    reverse, reverse!, copy, ==, inv, ≈, +, -, /
import DSP: conv

import InfiniteArrays: OrientedInfinity, Infinity

import PGFPlotsX: Plot, Options, PlotInc

export Infinity, ∞, InfiniteVector, CompactInfiniteVector, FixedInfiniteVector,  PeriodicInfiniteVector, downsample, upsample, δ, shift, shift!,
    ztransform, moment, fouriertransform, ⋆, ⊛, leastsquares_inv, eachnonzeroindex,
    alternating, alternating_flip, evenpart, oddpart, subvector

include("Integers.jl")





abstract type DoubleInfiniteArray{T,N} <: AbstractArray{T,N} end
Base.size(::DoubleInfiniteArray{T,N}) where {T,N} = ntuple(k->∞∞,Val(N))
Base.length(::DoubleInfiniteArray) = ∞∞
Base.eachindex(A::DoubleInfiniteArray) = CartesianIndices(axes(A))
Base.axes1(::DoubleInfiniteArray) = Integers()
Base.axes(::DoubleInfiniteArray{T,N}) where {T,N} = ntuple(k->Integers(),Val(N))
Base.checkbounds(::Type{Bool}, A::DoubleInfiniteArray, I...) = all(in.(Ref(eachindex(A)), I))

struct InfiniteArrayStyle{T,N} <: Broadcast.AbstractArrayStyle{N} end
Base.BroadcastStyle(::Type{DIA}) where {T,N,DIA<:DoubleInfiniteArray{T,N}} = InfiniteArrayStyle{DIA,N}()
Base.@pure function Broadcast.BroadcastStyle(a::InfiniteArrayStyle{DIA1,N}, b::InfiniteArrayStyle{DIA2,M}) where {N,M,S,T,DIA1<:DoubleInfiniteArray{S,N},DIA2<: DoubleInfiniteArray{T,M}}
    if Base.typename(DIA1) === Base.typename(DIA2)
        if M > N
            return b
        else
            return a
        end
    end
    return Broadcast.Unknown()
end
# Base.print_array(io::IO, X::DoubleInfiniteArray{T,2}) where T = Base.show_vector(io, X)
Base.LinearIndices(::NTuple{N,Integers}) where N = Integers()


const DoubleInfiniteVector{T} = DoubleInfiniteArray{T,1}
const InfiniteVector = DoubleInfiniteVector
include("DoubleInfiniteVector.jl")


include("AbstractCompactInfiniteVector.jl")
include("AbstractPeriodicInfiniteVector.jl")

include("CompactInfiniteVector.jl")


include("arithmetics.jl")
include("plots.jl")
# include("extensionsequences.jl")
# include("derivedsequences.jl")
# include("embeddingsequences.jl")

end # module
