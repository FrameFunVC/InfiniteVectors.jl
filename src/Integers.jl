struct Integers <: OrdinalRange{Int,Int}
end
(:)(a::OrientedInfinity{Bool}, b::Infinity) = a.angle == true ? Integers() : ∞:∞

Base.in(::Integers, ::Int) = true
Base.in(::Integers, ::Colon) = true
Base.in(::Integers, a) = false
Base.show(io::IO, ::Integers) = print(io, Symbol("\u2124"))

Base.length(::Integers) = ∞
Base.iterate(::Integers, a=nothing) = throw(MethodError("Can not iterate over infinite range"))
