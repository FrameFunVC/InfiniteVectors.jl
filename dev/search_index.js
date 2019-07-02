var documenterSearchIndex = {"docs":
[{"location":"#InfiniteVectors.jl-Documentation-1","page":"Home","title":"InfiniteVectors.jl Documentation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"A Julia package for infinite vectors operations such as convolution, z-transform.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"For installation instructions, see Installation.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"For a  full description of the functionality use the manual:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\"man/indices.md\",\n            \"man/doublyinfinite.md\"]","category":"page"},{"location":"#Installation-1","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"InfiniteVectors.jl is not added to General.jl. The package can easily be installed by cloning its git repository. From the Julia REPL, type ] to enter Pkg mode and run","category":"page"},{"location":"#","page":"Home","title":"Home","text":"pkg> add https://github.com/vincentcp/InfiniteVectors.jl","category":"page"},{"location":"#","page":"Home","title":"Home","text":"or in a file you could use","category":"page"},{"location":"#","page":"Home","title":"Home","text":"using Pkg\npkg\"add https://github.com/vincentcp/InfiniteVectors.jl\"","category":"page"},{"location":"#Development-1","page":"Home","title":"Development","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"For development instructions see the Development","category":"page"},{"location":"man/indices/#","page":"Indices","title":"Indices","text":"DocTestSetup = quote\n    using InfiniteVectors\nend\n","category":"page"},{"location":"man/indices/#Indices-1","page":"Indices","title":"Indices","text":"","category":"section"},{"location":"man/indices/#","page":"Indices","title":"Indices","text":"Even though an InfiniteVector is an AbstractVector","category":"page"},{"location":"man/indices/#","page":"Indices","title":"Indices","text":"julia> InfiniteVector <: AbstractVector\ntrue","category":"page"},{"location":"man/indices/#","page":"Indices","title":"Indices","text":"it does not support indexing since","category":"page"},{"location":"man/indices/#","page":"Indices","title":"Indices","text":"julia> vec = CompactInfiniteVector(1:3)\nCompactInfiniteVector{Int64} with indices ℤ:\n[  …, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, …  ]\njulia> I = eachindex(vec)\nℤ","category":"page"},{"location":"man/indices/#","page":"Indices","title":"Indices","text":"which is an object with infinite length.","category":"page"},{"location":"man/indices/#","page":"Indices","title":"Indices","text":"julia> length(I)\n∞","category":"page"},{"location":"man/indices/#","page":"Indices","title":"Indices","text":"note: Note\nThe function eachnonzeroindex is supported for vectors that have a finite number of non-zero elements.julia> eachnonzeroindex(vec)\n0:2","category":"page"},{"location":"man/vectors/#Vector-types-1","page":"Vector types","title":"Vector types","text":"","category":"section"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"This sections contains a list of the different supported vector types accompanied with the functions specific to those types.","category":"page"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"Pages = [\"vectors.md\"]\nDepth = 2","category":"page"},{"location":"man/vectors/#Compact-vectors-1","page":"Vector types","title":"Compact vectors","text":"","category":"section"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"The abstract type for vectors with a finite number of non-zero elements is","category":"page"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"InfiniteVectors.AbstractFiniteNZInfiniteVector","category":"page"},{"location":"man/vectors/#InfiniteVectors.AbstractFiniteNZInfiniteVector","page":"Vector types","title":"InfiniteVectors.AbstractFiniteNZInfiniteVector","text":"AbstractFiniteNZInfiniteVector{T} <: BiInfiniteVector{T}\n\nInstance of BiInfiniteVector that implements eachnonzeroindex(vec).\n\nAlso hascompactsupport(vec) == true\n\n\n\n\n\n","category":"type"},{"location":"man/vectors/#CompactInfiniteVector-1","page":"Vector types","title":"CompactInfiniteVector","text":"","category":"section"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"There are two types that support compact vectors. The second one is fully typed.","category":"page"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"CompactInfiniteVector\nFixedInfiniteVector","category":"page"},{"location":"man/vectors/#InfiniteVectors.CompactInfiniteVector","page":"Vector types","title":"InfiniteVectors.CompactInfiniteVector","text":"mutable struct CompactInfiniteVector{T} <: AbstractFiniteNZInfiniteVector{T}\n\nA CompactInfiniteVector contains a sequence of nonzero elements starting at a given offset.\n\nCompactInfiniteVector(a::AbstractVector{T}, offset = 0)\n\nConstruct a CompactInfiniteVector with indices starting at offset.\n\nExamples\n\njulia> CompactInfiniteVector(1:3,-1)\nCompactInfiniteVector{Int64} with indices ℤ:\n[  …, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, …  ]\n\n\n\n\n\n","category":"type"},{"location":"man/vectors/#InfiniteVectors.FixedInfiniteVector","page":"Vector types","title":"InfiniteVectors.FixedInfiniteVector","text":"struct FixedInfiniteVector{L,OFS,T} <: AbstractFiniteNZInfiniteVector{T}\n\nA FixedInfiniteVector is a fully typed sequence of nonzero elements starting at a given offset.\n\nFixedInfiniteVector(a::AbstractVector{T}, offset = 0)\n\nConstruct a FixedInfiniteVector with indices starting at offset.\n\nExamples\n\njulia> FixedInfiniteVector(1:3,-1)\nFixedInfiniteVector{3,-1,Int64} with indices ℤ:\n[  …, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, …  ]\n\n\n\n\n\n","category":"type"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"For the dirac delta a special constructor is present. ","category":"page"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"δ","category":"page"},{"location":"man/vectors/#InfiniteVectors.δ","page":"Vector types","title":"InfiniteVectors.δ","text":"δ(n::Int)\n\nConstruct Dirac delts, i.e., a bi-infinite vector with δ(k)=1, if k=n, and δ(k)=0, otherwise.\n\n\n\n\n\n","category":"function"},{"location":"man/vectors/#Functions-specific-to-compact-vectors-1","page":"Vector types","title":"Functions specific to compact vectors","text":"","category":"section"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"eachnonzeroindex\nsubvector\nsupport","category":"page"},{"location":"man/vectors/#InfiniteVectors.eachnonzeroindex","page":"Vector types","title":"InfiniteVectors.eachnonzeroindex","text":"eachnonzeroindex(vec)\n\nIndices of vec that contain non-zero elements. This function can be called if\n\nhascompactsupport(vec) == true\n\n\n\n\n\n","category":"function"},{"location":"man/vectors/#InfiniteVectors.subvector","page":"Vector types","title":"InfiniteVectors.subvector","text":"subvector(vec::CompactInfiniteVector) = vec.subvec\n\nThe vector of values at eachnonzeroindex\n\n\n\n\n\nsubvector(vec::AbstractPeriodicInfiniteVector)\n\nThe vector that goes from 0 to P-1, where P is the period of vec\"\n\n\n\n\n\n","category":"function"},{"location":"man/vectors/#InfiniteVectors.support","page":"Vector types","title":"InfiniteVectors.support","text":"support(vec::CompactInfiniteVector)\n\nThe minimum and maximum index of the non-zero elements\n\n\n\n\n\n","category":"function"},{"location":"man/vectors/#Periodic-vectors-1","page":"Vector types","title":"Periodic vectors","text":"","category":"section"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"InfiniteVectors.AbstractPeriodicInfiniteVector","category":"page"},{"location":"man/vectors/#InfiniteVectors.AbstractPeriodicInfiniteVector","page":"Vector types","title":"InfiniteVectors.AbstractPeriodicInfiniteVector","text":"abstract type AbstractPeriodicInfiniteVector{T} <: BiInfiniteVector{T} end\n\nA doubly infinite vector that has a period.\n\n\n\n\n\n","category":"type"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"PeriodicInfiniteVector","category":"page"},{"location":"man/vectors/#InfiniteVectors.PeriodicInfiniteVector","page":"Vector types","title":"InfiniteVectors.PeriodicInfiniteVector","text":"mutable struct PeriodicInfiniteVector{T} <: AbstractPeriodicInfiniteVector{T}\n\nA doubly infinite vector that has a period.\n\nPeriodicInfiniteVector(a::AbstractVector{T})\n\nConstruct a PeriodicInfiniteVector. The first element of a becomes the element at index zero, the second one at index 1.\n\nExamples\n\njulia> PeriodicInfiniteVector(1:3)\nPeriodicInfiniteVector{Int64} with indices ℤ:\n[  …, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, …  ]\n\n\n\n\n\n","category":"type"},{"location":"man/vectors/#Functions-specific-to-periodic-vectors-1","page":"Vector types","title":"Functions specific to periodic vectors","text":"","category":"section"},{"location":"man/vectors/#","page":"Vector types","title":"Vector types","text":"period","category":"page"},{"location":"man/vectors/#InfiniteVectors.period","page":"Vector types","title":"InfiniteVectors.period","text":"period(vec::AbstractPeriodicInfiniteVector)\n\nThe period of the periodic vector.\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#Vector-operations-1","page":"Vector operations","title":"Vector operations","text":"","category":"section"},{"location":"man/vectoroperations/#","page":"Vector operations","title":"Vector operations","text":"downsample\nupsample\nshift\nshift!\nreverse\nreverse!\nalternating_flip\nalternating\nevenpart\noddpart\ninv\nadjoint\ntranspose\nztransform\nfouriertransform\nmoment\nleastsquares_inv\n*\n⊛\n⋆","category":"page"},{"location":"man/vectoroperations/#InfiniteVectors.downsample","page":"Vector operations","title":"InfiniteVectors.downsample","text":"downsample(vec::InfiniteArrays, m::Int)\n\nThe resulting vector r satisfies r(k) = vec(mk)\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.upsample","page":"Vector operations","title":"InfiniteVectors.upsample","text":"upample(vec::InfiniteVector, m::Int)\n\nThe resulting vector r satisfies r(k) = vec(k/m), if k is multple of m, otherwise, r(k)=0\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.shift","page":"Vector operations","title":"InfiniteVectors.shift","text":"shift(vec::InfiniteArrays, m::Int)\n\nThe resulting vector r satisfies r(k) = vec(k+m)\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.shift!","page":"Vector operations","title":"InfiniteVectors.shift!","text":"shift!(vec::InfiniteArrays, m::Int)\n\nIn-place shifting of vector. (not always possible)\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#Base.reverse","page":"Vector operations","title":"Base.reverse","text":"reverse(vec::BiInfiniteVector)\n\nTime-reversel: vec(-k)\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#Base.reverse!","page":"Vector operations","title":"Base.reverse!","text":"reverse!(vec::CompactInfiniteVector)\n\nIn-place time reversel: vec(-k)\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.alternating_flip","page":"Vector operations","title":"InfiniteVectors.alternating_flip","text":"alternating_flip(vec::BiInfiniteVector, pivot = 1)\n\nFrom a given filter h(i), compute a new filter satisfying the alternating flip relation, centered around the given pivot:\n\ng(k) = (-1)^k h(pivot-k)\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.alternating","page":"Vector operations","title":"InfiniteVectors.alternating","text":"alternating(vec::InfiniteVector)\n\nCompute a new 'alternating' filter satisfying g(k) = (-1)^k h(k)\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.evenpart","page":"Vector operations","title":"InfiniteVectors.evenpart","text":"evenpart(vec::InfiniteVector)\n\nReturn the even part of a sequence `s`, defined by `s_e[k] = s[2k]`.\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.oddpart","page":"Vector operations","title":"InfiniteVectors.oddpart","text":"oddpart(vec::InfiniteVector)\n\nReturn the odd part of a sequence s, defined by s_o[k] = s[2k+1]\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#Base.inv","page":"Vector operations","title":"Base.inv","text":"inv(a::CompactInfiniteVector{T}, m::Int)\n\nA solution of a*b_m=δ, given a.\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#Base.adjoint","page":"Vector operations","title":"Base.adjoint","text":"adjoint(vec::BiInfiniteVector)\n\nReturns the paraconjugate overlinevec(-k)`\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#Base.transpose","page":"Vector operations","title":"Base.transpose","text":"transpose(vec::BiInfiniteVector)\n\nReturns the time-reversed vector vec(-k)\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.ztransform","page":"Vector operations","title":"InfiniteVectors.ztransform","text":"ztransform(vec::InfiniteVector, z)\n\nThe Z transform of a sequence is a continuous function of z, defined by S(z) = sum_kin ℤ s_k z^-k\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.fouriertransform","page":"Vector operations","title":"InfiniteVectors.fouriertransform","text":"fouriertransform(s::BiInfiniteVector, ω)\n\nThe Fourier transform of a sequence is defined by S(ω) = sum_kin ℤ s_k e^-i ω k. It is like the Z transform with z = e^i ω. The Fourier transform is a 2π-periodic continuous function of ω.\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.moment","page":"Vector operations","title":"InfiniteVectors.moment","text":"moment(vec::InfiniteVector, j)\n\nThe j-th discrete moment of a sequence is defined as sum_kin ℤ h_k k^j.\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.leastsquares_inv","page":"Vector operations","title":"InfiniteVectors.leastsquares_inv","text":"leastsquares_inv(v::BiInfiniteVector, m::Int)\n\nThe inverse filter f = v*(v*v_m)^-1_m\n\nreturns by applying on b v*f*b_m_m the least squares approximation of b in the range of v.\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#Base.:*","page":"Vector operations","title":"Base.:*","text":"*(a::InfiniteVector, b::InfiniteVector)\n\nThe convolution is defined as c(n) = sum_kin ℤ a(k)b(n-k)\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.:⊛","page":"Vector operations","title":"InfiniteVectors.:⊛","text":"⊛(a::AbstractPeriodicInfiniteVector, b::AbstractPeriodicInfiniteVector)\n\nThe circular convoluation defined as c(n) = sum_k=1^N a(k)b(n-k) for nin ℤ, where N is the period of a and b.\n\n\n\n\n\n","category":"function"},{"location":"man/vectoroperations/#InfiniteVectors.:⋆","page":"Vector operations","title":"InfiniteVectors.:⋆","text":"⋆(a::InfiniteVector, b::InfiniteVector)\n\nThe convolution is defined as c(n) = sum_kin ℤ a(k)b(n-k)\n\n\n\n\n\n","category":"function"},{"location":"development/#development-1","page":"Developer information","title":"Developer information","text":"","category":"section"},{"location":"development/#","page":"Developer information","title":"Developer information","text":"using InfiniteVectors","category":"page"},{"location":"development/#","page":"Developer information","title":"Developer information","text":"InfiniteVectors.ZTransform\nInfiniteVectors.Integers\nInfiniteVectors.BiInfiniteVector\nInfiniteVectors.hascompactsupport\nInfiniteVectors.nexteven\nInfiniteVectors.nextodd\nInfiniteVectors.previouseven\nInfiniteVectors.sublength","category":"page"},{"location":"development/#InfiniteVectors.ZTransform","page":"Developer information","title":"InfiniteVectors.ZTransform","text":"ZTransform{S<:BiInfiniteVector} <: Function\n\nZTransform is a wrapper type for the ztransform function.\n\n\n\n\n\n","category":"type"},{"location":"development/#InfiniteVectors.Integers","page":"Developer information","title":"InfiniteVectors.Integers","text":"struct InfiniteVectors.Integers <: OrdinalRange{Int,Int}\n\nStruct representing all integers.\n\n\n\n\n\n","category":"type"},{"location":"development/#InfiniteVectors.BiInfiniteVector","page":"Developer information","title":"InfiniteVectors.BiInfiniteVector","text":"Any subtype of InfiniteVector is a bi-infinite list of values, with one value for each integer in Z.\n\nThey are not iterable\n\n\n\n\n\n","category":"type"},{"location":"development/#InfiniteVectors.hascompactsupport","page":"Developer information","title":"InfiniteVectors.hascompactsupport","text":"hascompactsupport(vec::InfiniteVector)\n\nAre their a finite number of non-zero elements.\n\n\n\n\n\n","category":"function"},{"location":"development/#InfiniteVectors.nexteven","page":"Developer information","title":"InfiniteVectors.nexteven","text":"The first even number greater than or equal to n.\n\n\n\n\n\n","category":"function"},{"location":"development/#InfiniteVectors.nextodd","page":"Developer information","title":"InfiniteVectors.nextodd","text":"The first odd number greater than or equal to n.\n\n\n\n\n\n","category":"function"},{"location":"development/#InfiniteVectors.previouseven","page":"Developer information","title":"InfiniteVectors.previouseven","text":"The last even number, smaller than or equal to n.\n\n\n\n\n\n","category":"function"},{"location":"development/#InfiniteVectors.sublength","page":"Developer information","title":"InfiniteVectors.sublength","text":"sublength(vec::InfiniteVector) = vec.subvec\n\nThe number of non-zero indices\n\n\n\n\n\n","category":"function"}]
}