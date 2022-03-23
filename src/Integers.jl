"""
   DoubleInfinity()
represents double infinite cardinality. Note that `DoubleInfinity <: Integer` to support
being treated as an index.
"""
struct DoubleInfinity <: Integer end
export ∞∞
const ∞∞ = DoubleInfinity()

Base.show(io::IO, ::DoubleInfinity) = print(io, "∞∞")
Base.string(::DoubleInfinity) = "∞∞"
Base.isinf(::DoubleInfinity) = true
Base.isfinite(::DoubleInfinity) = false

import Base: ==, *, -, +, fld, cld, div,≤,<,>,≥
import InfiniteArrays: NotANumber
==(x::DoubleInfinity, y::DoubleInfinity) = true
Base.sign(y::DoubleInfinity) = 1
Base.angle(x::DoubleInfinity) = 0
Base.zero(::DoubleInfinity) = 0

==(x::DoubleInfinity, y::Number) = isinf(y) && angle(y) == angle(x)
==(y::Number, x::DoubleInfinity) = x == y

Base.isless(x::DoubleInfinity, y::DoubleInfinity) = false
Base.isless(x::Real, y::DoubleInfinity) = isfinite(x) || sign(y) == -1
Base.isless(x::AbstractFloat, y::DoubleInfinity) = isless(x, convert(typeof(x), y))
Base.isless(x::DoubleInfinity, y::AbstractFloat) = false
Base.isless(x::DoubleInfinity, y::Real) = false

+(::DoubleInfinity, ::DoubleInfinity) = ∞∞
+(::Number, y::DoubleInfinity) = ∞∞
+(::DoubleInfinity, ::Number) = ∞∞
-(::DoubleInfinity, ::Number) = ∞∞
-(x::Number, ::DoubleInfinity) = x + (-∞∞)

+(::Integer, y::DoubleInfinity) = ∞∞
+(::DoubleInfinity, ::Integer) = ∞∞
-(::DoubleInfinity, ::Integer) = ∞∞
-(x::Integer, ::DoubleInfinity) = x + (-∞∞)
+(::Complex, y::DoubleInfinity) = ∞∞
+(::DoubleInfinity, ::Complex) = ∞∞
-(::DoubleInfinity, ::Complex) = ∞∞
-(x::Complex, ::DoubleInfinity) = x + (-∞∞)

-(::DoubleInfinity, ::DoubleInfinity) = NotANumber()

# ⊻ is xor
*(::DoubleInfinity) = ∞∞
*(::DoubleInfinity, ::DoubleInfinity) = ∞∞



for OP in (:fld,:cld,:div)
  @eval begin
    $OP(::DoubleInfinity, ::Real) = fld
    $OP(::DoubleInfinity, ::DoubleInfinity) = NotANumber()
  end
end

div(::T, ::DoubleInfinity) where T<:Real = zero(T)
fld(x::T, ::DoubleInfinity) where T<:Real = signbit(x) ? -one(T) : zero(T)
cld(x::T, ::DoubleInfinity) where T<:Real = signbit(x) ? zero(T) : one(T)

Base.mod(::DoubleInfinity, ::DoubleInfinity) = NotANumber()
Base.mod(::DoubleInfinity, ::Real) = NotANumber()
function Base.mod(x::Real, ::DoubleInfinity)
    x ≥ 0 || throw(ArgumentError("mod(x,∞) is unbounded for x < 0"))
    x
end

Base.min(::DoubleInfinity, ::DoubleInfinity) = ∞
Base.max(::DoubleInfinity, ::DoubleInfinity) = ∞
Base.min(x::Real, ::DoubleInfinity) = x
Base.max(::Real, ::DoubleInfinity) = ∞
Base.min(::DoubleInfinity, x::Real) = x
Base.max(::DoubleInfinity, ::Real) = ∞

≤(::DoubleInfinity, ::DoubleInfinity) = true
<(::DoubleInfinity, ::DoubleInfinity) = false
≥(::DoubleInfinity, ::DoubleInfinity) = true
>(::DoubleInfinity, ::DoubleInfinity) = false

for OP in (:<, :≤)
    @eval begin
        $OP(::Real, ::DoubleInfinity) = true
        $OP(::DoubleInfinity, ::Real) = false
    end
end

for OP in (:>, :≥)
    @eval begin
        $OP(::Real, ::DoubleInfinity) = false
        $OP(::DoubleInfinity, ::Real) = true
    end
end

"""
    struct InfiniteVectors.Integers{I} <: OrdinalRange{I,I}

Struct representing all integers.
"""
struct Integers{I} <: AbstractUnitRange{I}
end
Integers() = Integers{Int}()

Base.in(::Integers, ::Int) = true
Base.in(::Integers, ::Colon) = true
Base.in(::Integers, a) = false
Base.show(io::IO, ::Integers) = print(io, Symbol("\u2124"))

Base.length(::Integers) = ∞∞
Base.iterate(::Integers, a=nothing) = throw(MethodError("Can not iterate over infinite range"))
Base.first(::Integers) = -∞
Base.last(::Integers) = +∞
Base.axes1(i::Integers) = i
Base.axes(i::Integers) = (i,)
