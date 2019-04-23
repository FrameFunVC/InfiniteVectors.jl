# sequence.jl

__precompile__()
module Sequences
using InfiniteArrays, DSP, LinearAlgebra, FFTW, PGFPlotsX

using StaticArrays

import Base: size, length, getindex, setindex!, similar, adjoint, transpose, sum, *, eltype, :,
    reverse, reverse!, copy, ==, inv, ≈
import DSP: conv

import InfiniteArrays: OrientedInfinity, Infinity

import PGFPlotsX: Plot, Options

# export evenpart, oddpart, alternating_flip, alternating
#
# export firstindex, lastindex, each_nonzero_index

# # Extension sequences
# export PeriodicExtension, ZeroPadding, ConstantPadding, ShiftedExtension,
#     SymmetricExtension, UndefinedExtension
#
# export each_subindex, subvector, sublength
#
# # Derived sequences
# export UpsampledInfiniteVector, DownsampledInfiniteVector, ReversedInfiniteVector, ShiftedInfiniteVector
#
# # Embedding sequences
# export EmbeddingInfiniteVector, PeriodicEmbedding, SymmetricEmbedding, FunctionEmbedding

export Infinity, ∞, CompactInfiniteVector, PeriodicInfiniteVector, downsample, upsample, δ, shift, shift!,
    ztransform, moment, fouriertransform, *, ⋆, ⊛, hascompactsupport, period, inv

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
