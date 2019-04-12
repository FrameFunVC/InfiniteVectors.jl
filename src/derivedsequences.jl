# derivedsequences.jl


######################################
# Type hierarchy and general interface
######################################

"""
A DerivedInfiniteVector is any sequence that is derived from another sequence.
Examples include a DownsampledInfiniteVector and an UpsampledInfiniteVector.
DerivedInfiniteVectors are lazy. No computation is performed until the sequence is
evaluated. This can be done for compactly supported sequences using `collect`.
"""
abstract type DerivedInfiniteVector{S} <: InfiniteVector end

eltype{S}(::Type{DerivedInfiniteVector{S}}) = eltype(S)
eltype{DS <: DerivedInfiniteVector}(::Type{DS}) = eltype(super(DS))

# If the subtypes only transform indices, then they can implement mapindex and imapindex,
# where mapindex maps an index 'k' of the derived sequence into an index 'l' of the original sequence,
# and imapindex is the inverse.
# We always use the symbol 'l' for an index of the original sequence and 'k' for the derived sequence.

# Default indexing rules are:
getindex{M}(s::DerivedInfiniteVector{M}, k) = s.seq[mapindex(s, k)]

setindex!{M}(s::DerivedInfiniteVector{M}, val, k) = s.seq[mapindex(s, k)] = val

"Return the original sequence of the derived sequence."
sequence(s::DerivedInfiniteVector) = s.seq

# Traits
hascompactsupport{S}(::Type{DerivedInfiniteVector{S}}) = hascompactsupport(S)


collect(s::DerivedInfiniteVector) = _collect(s, hascompactsupport(s))

_collect(s::DerivedInfiniteVector, compactsupport::True) = CompactInfiniteVector(eltype(s)[s[i] for i in fistindex(s):lastindex(s)], firstindex(s))

_collect(s::DerivedInfiniteVector, compactsupport::False) = throw(BoundsError())



######################################
# Definitions of derived sequences
######################################


"""
A DownsampledInfiniteVector is determined by a sequence `s`, a downsampling factor `M`
and a `shift`. It is defined by `ds_k = s_{shift+M*k}`.

A DownsampledInfiniteVector acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 's' will be modified.
"""
struct DownsampledInfiniteVector{M,S} <: DerivedInfiniteVector{S}
    seq     ::  S
    shift   ::  Int

    DownsampledInfiniteVector{M,S}(seq::InfiniteVector, shift) where {M,S} = new(seq, shift)
end

# Default downsampling factor is 2.
DownsampledInfiniteVector(s) = DownsampledInfiniteVector(s, Val{2})

# This one is convenient but not type-stable.
DownsampledInfiniteVector(s, m::Int, shift) = DownsampledInfiniteVector(s, Val{m}, shift)

# Default shift is 0.
DownsampledInfiniteVector{M,S}(s::S, ::Type{Val{M}}, shift = 0) = DownsampledInfiniteVector{M,S}(s, shift)

downsample(s::InfiniteVector, args...) = DownsampledInfiniteVector(s, args...)


mapindex{M}(s::DownsampledInfiniteVector{M}, k) = s.shift+M*k

imapindex{M}(s::DownsampledInfiniteVector{M}, l) = mod(l-s.shift, M) == 0 ? div(l-s.shift, M) : throw(BoundsError())



"""
An UpsampledInfiniteVector is determined by a sequence `s`, an upsampling factor `M`
and a `shift`. It is defined by `us_k = s_l`, for `k = shift+M*l`, and `us_k = 0` otherwise.

An UpsampledInfiniteVector acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 's' will be modified.
"""
struct UpsampledInfiniteVector{M,S} <: DerivedInfiniteVector{S}
    seq     ::  S
    shift   ::  Int

    UpsampledInfiniteVector{M,S}(seq::InfiniteVector, shift) where {M,S} = new(seq, shift)
end

# Default upsampling factor is 2.
UpsampledInfiniteVector(s) = UpsampledInfiniteVector(s, Val{2})

# Default shift is 0.
UpsampledInfiniteVector{M,S}(s::S, ::Type{Val{M}}, shift = 0) = UpsampledInfiniteVector{M,S}(s, shift)

# This one is convenient but not type-stable.
UpsampledInfiniteVector(s, m::Int, shift) = UpsampledInfiniteVector(s, Val{m}, shift)

upsample(s::InfiniteVector, args...) = UpsampledInfiniteVector(s, args...)


mapindex{M}(s::UpsampledInfiniteVector{M}, k) = mod(k-s.shift, M) == 0 ? k : throw(BoundsError())

getindex{M}(s::UpsampledInfiniteVector{M}, k) = mod(k-s.shift, M) == 0 ? s.seq[k] : zero(eltype(s))

each_nonzero_index{M}(s::UpsampledInfiniteVector{M}) =
    s.shift+imapindex(s, firstindex(s.seq)):M:s.shift+imapindex(s, lastindex(s.seq))


"""
A ReversedInfiniteVector is the time-reversal of a given sequence `s`, i.e. `rs_k = s_{-k}`.

A ReversedInfiniteVector acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 's' will be modified.
"""
struct ReversedInfiniteVector{S} <: DerivedInfiniteVector{S}
    seq     ::  S
end

mapindex(s::ReversedInfiniteVector, k) = -k

imapindex(s::ReversedInfiniteVector, l) = -l

firstindex(s::ReversedInfiniteVector) = -lastindex(sequence(s))

lastindex(s::ReversedInfiniteVector) = -firstindex(sequence(s))

"Reverse the given sequence."
reverse(s::InfiniteVector) = ReversedInfiniteVector(s)

reverse(s::ReversedInfiniteVector) = sequence(s)




"""
A ShiftedInfiniteVector is determined by a `shift` and a given sequence `s`. It is
defined by `ss_k = s_{k+shift}`.

A ShiftedInfiniteVector acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 's' will be modified.
"""
struct ShiftedInfiniteVector{S} <: DerivedInfiniteVector{S}
    seq     ::  S
    shift   ::  Int
end

mapindex(s::ShiftedInfiniteVector, k) = k + s.shift

impaindex(s::ShiftedInfiniteVector, l) = l - s.shift

"Shift a sequence forward by `shift` positions."
shift(s::InfiniteVector, shift::Int) = ShiftedInfiniteVector(s, shift)

shift(s::ShiftedInfiniteVector, shift::Int) = ShiftedInfiniteVector(sequence(s), shift+s.shift)

firstindex(s::ShiftedInfiniteVector) = firstindex(sequence(s)) + s.shift

lastindex(s::ShiftedInfiniteVector) = lastindex(sequence(s)) + s.shift
