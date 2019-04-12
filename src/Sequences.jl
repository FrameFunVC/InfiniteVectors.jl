# sequence.jl

__precompile__()
module Sequences
using InfiniteArrays, DSP, LinearAlgebra
using DomainSets: GeometricSpace

using StaticArrays

import Base: size, length, getindex, setindex!, similar, adjoint, transpose, sum, *, eltype,
    reverse
import DSP: conv

#, eachindex, collect, &, |, *, transpose, ctranspose, conj, sum, +, -, /,
#             convert, widen, reverse

# # Main abstract types
# export InfiniteVector, ExtensionInfiniteVector, DerivedInfiniteVector
#
# # Utility function
# export promote_eltype
# # Traits
# export True, False, hascompactsupport
#
# export moment, ztransform, fouriertransform, support
#
# export evenpart, oddpart, alternating_flip, alternating
#
# export firstindex, lastindex, each_nonzero_index
#
# # Traits
# export True, False, hascompactsupport

# # Extension sequences
# export PeriodicExtension, ZeroPadding, ConstantPadding, ShiftedExtension,
#     SymmetricExtension, UndefinedExtension
#
# export each_subindex, subvector, sublength
#
# # Compactly supported sequences
# export AbstractCompactInfiniteVector, FixedInfiniteVector
#
# # Derived sequences
# export UpsampledInfiniteVector, DownsampledInfiniteVector, ReversedInfiniteVector, ShiftedInfiniteVector
#
# # Embedding sequences
# export EmbeddingInfiniteVector, PeriodicEmbedding, SymmetricEmbedding, FunctionEmbedding
#
# export shift, reverse, upsample, downsample
#
export Infinity, CompactInfiniteVector, PeriodicInfiniteVector


include("AbstractDoubleInfiniteVector.jl")
include("AbstractCompactInfiniteVector.jl")
include("AbstractPeriodicInfiniteVector.jl")

include("CompactInfiniteVector.jl")


# include("compactsequences.jl")
# include("extensionsequences.jl")
# include("derivedsequences.jl")
# include("embeddingsequences.jl")
#
#
# promote_eltype{T}(s::AbstractCompactInfiniteVector{T}, ::Type{T}) = s

end # module
