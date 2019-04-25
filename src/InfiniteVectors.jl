# sequence.jl

__precompile__()
module InfiniteVectors
using InfiniteArrays, DSP, LinearAlgebra, FFTW, PGFPlotsX, FastTransforms

using StaticArrays

import Base: size, length, getindex, setindex!, similar, adjoint, transpose, sum, *, eltype, :,
    reverse, reverse!, copy, ==, inv, ≈
import DSP: conv

import InfiniteArrays: OrientedInfinity, Infinity

import PGFPlotsX: Plot, Options, PlotInc

export Infinity, ∞, CompactInfiniteVector, PeriodicInfiniteVector, downsample, upsample, δ, shift, shift!,
    ztransform, moment, fouriertransform, *, ⋆, ⊛, hascompactsupport, period, inv, leastsquares_inv

include("Integers.jl")


abstract type InfiniteVector{T} <: AbstractVector{T} end
Base.size(::InfiniteVector) = (∞,)
Base.length(::InfiniteVector) = ∞



include("AbstractDoubleInfiniteVector.jl")
include("AbstractCompactInfiniteVector.jl")
include("AbstractPeriodicInfiniteVector.jl")

include("CompactInfiniteVector.jl")


include("arithmetics.jl")
include("plots.jl")
# include("extensionsequences.jl")
# include("derivedsequences.jl")
# include("embeddingsequences.jl")

end # module
