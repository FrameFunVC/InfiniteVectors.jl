using Sequences, LinearAlgebra, Test, DSP

struct MyConcreteInfiniteVector{T} <: Sequences.AbstractDoubleInfiniteVector{T}
end

f = MyConcreteInfiniteVector{Float64}();
@test_throws MethodError iterate(eachindex(f))
@test_throws MethodError  iterate(f)
@test size(f) == (∞,)
@test length(f) == ∞



function test_inf_vector(a)
    b = copy(a)
    reverse!(b)
    reverse!(b)
    @test a==b
    for m in 1:12
        @test a[-10m:10m][1:m:end]==downsample(a, m)[-10:10]
        @test a==downsample(upsample(a, m), m)
        @test reverse(a)[-10:10] == a[10:-1:-10]

        @test a' == reverse(a)
        @test (δ(3) * a) == shift(a,3)
        b = copy(a)
        @test shift(a, m) == shift!(b, m)
    end
end

@testset "basis functionality of CompactInfiniteVector" begin
    for offset in -4:4
        a = CompactInfiniteVector(1:9,offset)
        test_inf_vector(a)
    end
end

@testset "basic functionality of PeriodicInfiniteVector"
    for i in 1:10
        a = PeriodicInfiniteVector(1:i)
        test_inf_vector(a)
    end
end
