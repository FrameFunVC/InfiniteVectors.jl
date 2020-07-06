using InfiniteVectors, LinearAlgebra, Test, DSP, PGFPlotsX


struct MyConcreteInfiniteVector{T} <: InfiniteVectors.DoubleInfiniteVector{T}
end

@testset "generic tests" begin
    f = MyConcreteInfiniteVector{Float64}();
    Base.getindex(::MyConcreteInfiniteVector, i) = 1.
    @test_throws MethodError iterate(eachindex(f))
    @test_throws MethodError  iterate(f)
    @test size(f) == (∞∞,)
    @test length(f) == ∞∞
    @test axes(f)===tuple(Base.axes1(f))
    io = IOBuffer()
    show(io,f)
    @test String(take!(io))=="[  …, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, …  ]"
end

function test_inf_vector(a; inplace=true)
    b = copy(a)
    if inplace
        reverse!(b)
        reverse!(b)
    end
    @test a==b
    for m in 1:12
        @test a[-10m:10m][1:m:end]==downsample(a, m)[-10:10]
        @test a==downsample(upsample(a, m), m)
        @test reverse(a)[-10:10] == a[10:-1:-10]
        @test alternating(reverse(a))[-10:10] == alternating_flip(a,0)[-10:10]
        @test evenpart(a)[-4:4] == a[-8:2:8]
        @test oddpart(a)[-4:4] == a[-7:2:9]

        @test a' == reverse(a)
        @test transpose(a) == reverse(a)
        @test (δ(3) * a) == shift(a,3)
        @test (a ⋆ δ(3)) == shift(a,3)
        b = copy(a)
        if inplace
            @test shift(a, m) == shift!(b, m)
        end
    end
end

@testset "basis functionality of CompactInfiniteVector" begin
    for offset in -4:4
        a = CompactInfiniteVector(1:9,offset)
        test_inf_vector(a)
    end

    # ztransform
    a = CompactInfiniteVector(ones(3),-1)
    @test sum(a) ≈ 3
    omega = LinRange(-1,1,10).*2pi
    x = exp.(1im*omega)
    @test ztransform.(Ref(upsample(a,2)),x) ≈ 1 .+ 2cos.(2omega)
    @test fouriertransform.(Ref(upsample(a,2)),omega) ≈ 1 .+ 2cos.(2omega)
    @test moment(a, 2) ≈ 2
    @test moment(upsample(a,2), 2) ≈ 8
    @test InfiniteVectors.hascompactsupport(a)
    c = PeriodicInfiniteVector(CompactInfiniteVector(1:4,-3),2)
    @test period(c) == 2
    @test c[0:1] == [6,4]

end

@testset "basic functionality of PeriodicInfiniteVector" begin
    for i in 1:10
        a = PeriodicInfiniteVector(1:i)
        test_inf_vector(a; inplace=false)
    end
    a = PeriodicInfiniteVector(1.:1.:5.)
    b = PeriodicInfiniteVector(ones(2))
    @test sum(InfiniteVectors.subvector(a)) ≈ sum(1:5)
    @test !InfiniteVectors.hascompactsupport(a)
    @test a⊛b isa PeriodicInfiniteVector
end

@testset "basis functionality of FixedInfiniteVector" begin
    for offset in -4:4
        a = FixedInfiniteVector(1:9,offset)
        test_inf_vector(a; inplace=false)
    end

    # ztransform
    a = FixedInfiniteVector(ones(3),-1)
    @test sum(a) ≈ 3
    omega = LinRange(-1,1,10).*2pi
    x = exp.(1im*omega)
    @test ztransform.(Ref(upsample(a,2)),x) ≈ 1 .+ 2cos.(2omega)
    @test fouriertransform.(Ref(upsample(a,2)),omega) ≈ 1 .+ 2cos.(2omega)
    @test moment(a, 2) ≈ 2
    @test moment(upsample(a,2), 2) ≈ 8
    @test InfiniteVectors.hascompactsupport(a)
end

@testset "basis functionality of CompactPeriodicInfiniteVector" begin
    for offset in -4:4, p in [7,8,13,14]
        a = CompactPeriodicInfiniteVector(1:9,p,offset)
        test_inf_vector(a; inplace=false)
    end
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
        @test downsample(a*inv(a, m;K=n>>1), m)[-10:10] ≈ δ(0)[-10:10]
    end
    io = IOBuffer()
    show(io, a)
end
