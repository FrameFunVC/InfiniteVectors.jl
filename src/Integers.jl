"""
    struct InfiniteVectors.Integers{I} <: OrdinalRange{I,I}

Struct representing all integers.
"""
struct Integers{I} <: OrdinalRange{I,I}
end
Integers() = Integers{Int}()
(:)(a::OrientedInfinity{Bool}, b::Infinity) = a.angle == true ? Integers() : ∞:∞

Base.in(::Integers, ::Int) = true
Base.in(::Integers, ::Colon) = true
Base.in(::Integers, a) = false
Base.show(io::IO, ::Integers) = print(io, Symbol("\u2124"))

Base.length(::Integers) = ∞
Base.iterate(::Integers, a=nothing) = throw(MethodError("Can not iterate over infinite range"))
Base.first(::Integers{I}) where I = typemin(I)
Base.last(::Integers{I}) where I = typemax(I)
