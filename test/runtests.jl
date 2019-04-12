using Sequences, LinearAlgebra, Test

struct MyConcreteInfiniteVector{T} <: Sequences.AbstractDoubleInfiniteVector{T}
end

f = MyConcreteInfiniteVector{Float64}();
@test_throws MethodError iterate(eachindex(f))
@test_trhows MethodError = iterate(f)
@test size(f) == (∞,)
@test length(f) == ∞
