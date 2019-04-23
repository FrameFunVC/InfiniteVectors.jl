using Sequences, LinearAlgebra, Test, DSP, PGFPlotsX

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

    # ztransform
    a = CompactInfiniteVector(ones(3),-1)
    omega = LinRange(-1,1,10).*2pi
    x = exp.(1im*omega)
    @test ztransform.(Ref(upsample(a,2)),x) ≈ 1 .+ 2cos.(2omega)
    @test fouriertransform.(Ref(upsample(a,2)),omega) ≈ 1 .+ 2cos.(2omega)
    @test moment(a, 2) ≈ 2
    @test moment(upsample(a,2), 2) ≈ 8
end

@testset "basic functionality of PeriodicInfiniteVector" begin
    for i in 1:10
        a = PeriodicInfiniteVector(1:i)
        test_inf_vector(a)
    end
    a = PeriodicInfiniteVector(1.:1.:5.)
    b = PeriodicInfiniteVector(ones(2))
    @test sum(Sequences.subvector(a)) ≈ sum(1:5)
end

@testset "plot" begin
    Plot(PeriodicInfiniteVector(rand(10)))
    Plot(PeriodicInfiniteVector(rand(ComplexF64, 10)))
    @pgf Plot({samples_at=-2:2}, PeriodicInfiniteVector(rand(10)))
    @pgf Plot({samples_at="1,2,...,10"}, PeriodicInfiniteVector(rand(10)))
end

@testset "inv" begin
    a = PeriodicInfiniteVector(rand(10))
    @test (a*inv(a))[0:9] ≈ δ(0)[0:9]
    a = PeriodicInfiniteVector(rand(ComplexF64, 10))
    @test (a*inv(a))[0:9] ≈ δ(0)[0:9]
    for n in 2:7, os in -3:3, m in 2:5
        a = CompactInfiniteVector(rand(n),os)
        @test downsample(a*inv(a, m), m)[-10:10] ≈ δ(0)[-10:10]
    end
end
