module GaussianRandomVariables

using ThickNumbers
using SpecialFunctions: SpecialFunctions, gamma

import Base: +, -, *, /, //, ^, inv
import Base: abs, abs2, max, min, sin, cos, sincos, sqrt
import Base: log, exp

export GVar, ±

const BaseReals = Union{AbstractFloat, Integer, AbstractIrrational, Rational}

## Type and constructors

struct GVar{T<:Real} <: ThickNumber{T}
    center::T
    σ::T

    function GVar{T}(center::Real, σ::Real) where T<:Real
        # σ < zero(σ) && throw(DomainError(σ, "σ must be nonnegative"))
        new{T}(center, σ)
    end
end

GVar{T}(x::Real) where T = GVar{T}(x, zero(x))
GVar(x::T) where T<:Real = GVar{T}(x)

GVar(center::T, σ::T) where T<:Real = GVar{T}(center, σ)
GVar(center::Real, σ::Real) = GVar(promote(center, σ)...)

GVar(center::T, σ::T) where T<:Integer = GVar(float(center), float(σ))
GVar(center::T, σ::T) where T<:Irrational = GVar(float(center), float(σ))

GVar(x::GVar) = x
GVar{T}(x::GVar) where T = GVar{T}(x.center, x.σ)

±(center::Real, σ::Real) = GVar(center, σ)

Base.promote_rule(::Type{GVar{T}}, ::Type{GVar{S}}) where {T<:Real,S<:Real} = GVar{promote_type(T,S)}
Base.promote_rule(::Type{GVar{T}}, ::Type{S}) where {T<:Real,S<:BaseReals} = GVar{promote_type(T,S)}

AbstractFloat(a::GVar{<:AbstractFloat}) = a
AbstractFloat(a::GVar) = GVar(AbstractFloat(a.center), AbstractFloat(a.σ))

## ThickNumbers API

ThickNumbers.loval(a::GVar) = a.center - a.σ
ThickNumbers.hival(a::GVar) = a.center + a.σ
ThickNumbers.lohi(::Type{G}, lo, hi) where G<:GVar = G((lo + hi)/2, nextfloat((hi - lo)/2))
ThickNumbers.midrad(::Type{G}, center, σ) where G<:GVar = G(center, σ)
ThickNumbers.basetype(::Type{GVar{T}}) where T = GVar
ThickNumbers.basetype(::Type{GVar}) = GVar
ThickNumbers.emptyset(::Type{G}) where G<:GVar = G(0, -1)

## Display

Base.show(io::IO, a::GVar) = print(io, a.center, " ± ", a.σ)


## Trait functions and constants

Base.zero(::GVar{T}) where T<:Real = zero(GVar{T})
Base.zero(::Type{GVar{T}}) where T<:Real = GVar(zero(T), zero(T))
Base.oneunit(::GVar{T}) where T<:Real = oneunit(GVar{T})
Base.oneunit(::Type{GVar{T}}) where T<:Real = GVar(oneunit(T),zero(T))

Base.real(a::GVar) = a
Base.conj(a::GVar) = a

function Base.hash(x::GVar, h::UInt)
    magic = Int === Int64 ? 0x21f1b1afbb07c31a : 0x237f8305
    h += magic
    return hash(x.center, hash(x.σ, h))
end


## Arithmetic

ispos(x::Real) = x > zero(x)
isneg(x::Real) = x < zero(x)
isnonneg(x::Real) = x >= zero(x)
isnonpos(x::Real) = x <= zero(x)

## Addition and subtraction
function +(a::GVar{<:Real}, b::Real)
    GVar(a.center + b, a.σ)
end
+(a::GVar{<:Integer}, b::Integer) = float(a) + float(b)
+(b::Real, a::GVar{<:Real}) = a+b

function -(a::GVar{<:Real}, b::Real)
    GVar(a.center - b, a.σ)
end
-(a::GVar{<:Integer}, b::Integer) = float(a) - float(b)
function -(b::Real, a::GVar{<:Real})
    GVar(b - a.center, a.σ)
end
-(b::Integer, a::GVar{<:Integer}) = float(b) - float(a)

function +(a::GVar{<:Real}, b::GVar{<:Real})
    @fastmath begin
        ret = GVar(a.center + b.center, sqrt(a.σ^2 + b.σ^2))
        return (isempty(a) | isempty(b)) ? emptyset(ret) : ret
    end
end
+(a::GVar{<:Integer}, b::GVar{<:Integer}) = float(a) + float(b)

function -(a::GVar{<:Real}, b::GVar{<:Real})
    @fastmath begin
        ret = GVar(a.center - b.center, sqrt(a.σ^2 + b.σ^2))
        return (isempty(a) | isempty(b)) ? emptyset(ret) : ret
    end
end
-(a::GVar{<:Integer}, b::GVar{<:Integer}) = float(a) - float(b)

## Multiplication
# Restrict to BaseReals so we avoid ambiguities with ForwardDiff
function *(x::BaseReals, a::GVar{<:Real})
    @fastmath begin
        c, σ = x*a.center, abs(x)*a.σ
        return GVar(c, σ)
    end
end
*(a::GVar{<:Real}, x::BaseReals) = x*a
# Prevent overflow by promoting to float
*(x::Integer, a::GVar{<:Integer}) = float(x)*float(a)
*(a::GVar{<:Integer}, x::Integer) = float(x)*float(a)

function *(a::GVar{<:Real}, b::GVar{<:Real})
    @fastmath begin
        ret = GVar(a.center*b.center, sqrt(a.center^2*b.σ^2 + b.center^2*a.σ^2 + a.σ^2*b.σ^2))
        return (isempty(a) | isempty(b)) ? emptyset(ret) : ret
    end
end
*(a::GVar{<:Integer}, b::GVar{<:Integer}) = float(a)*float(b)   # prevent overflow

## Division
function /(a::GVar{<:Real}, x::Real)
    @fastmath begin
        c, σ = a.center/x, a.σ / abs(x)
        return GVar(c, σ)
    end
end

function inv(a::GVar{T}) where T<:AbstractFloat  # must be AbstractFloat so typemax gives Inf
    z = zero(1/oneunit(T))
    isempty(a) && return emptyset(basetype(a){typeof(z)})
    zero(T) ∈ a && return GVar(iszero(a.center) ? z : 1/a.center, typemax(T))
    # As of 2023-10-11, the explicit case for `var[X/Y]` on
    #    https://en.wikipedia.org/wiki/Taylor_expansions_for_the_moments_of_functions_of_random_variables#Second_moment
    # is badly wrong. Use the general case below.
    return GVar(meanvar(inv, x -> -1/x^2, x -> 2/x^3, a.center, a.σ)...)
end
inv(a::GVar{<:Real}) = inv(float(a))


/(a::Real, b::GVar{<:Real}) = a*inv(b)

/(a::GVar{<:Real}, b::GVar{<:Real}) = a*inv(b)

//(a::GVar, b::GVar) = a / b    # to deal with rationals

## Powers
Base.literal_pow(::typeof(^), x::GVar, ::Val{0}) = oneunit(x)
Base.literal_pow(::typeof(^), x::GVar, ::Val{1}) = x
function Base.literal_pow(::typeof(^), x::GVar, ::Val{2})
    c2, σ2 = x.center^2, x.σ^2
    return GVar(c2 + σ2, sqrt(4*c2*σ2 + 2*σ2^2)*sign(x.σ))
end
Base.literal_pow(::typeof(^), x::GVar, ::Val{p}) where p = x^p
function ^(a::GVar, p::Integer)
    isempty(a) && return a
    p < 0 && return inv(a^(-p))
    cp2 = a.center ^ (p-2)
    c2, σ2 = a.center^2, a.σ^2
    # If |c| ≫ σ, then these are accurate:
    c′ = c2*cp2 + p*(p-1)*σ2*cp2/2
    σ′² = p^2*σ2*cp2^2*c2
    # But if |c| ≪ σ and p is even, the dominant term is from σ^p. So let's just add that in:
    if iseven(p)
        σ′² += (gamma(p + 0.5) * 2^p / sqrt(pi) - 1) * a.σ^(2p)  # the -1 is from subtracting <a^p>^2
    end
    return GVar(c′, sqrt(σ′²))
end

function ^(a::GVar, p::Real)
    isinteger(p) && return a^Int(p)
    isempty(a) && return a
    return GVar(meanvar(x->x^p, x->p*x^(p-1), x -> p*(p-1)*x^(p-2), a.center, a.σ)...)
end

## Functions
function meanvar(f::F, df::DF, ddf::DDF, c, σ) where {F, DF, DDF}
    # As of 2023-10-11, the Wikipedia page
    #   https://en.wikipedia.org/wiki/Taylor_expansions_for_the_moments_of_functions_of_random_variables#Second_moment
    # is quite confused about this: it presents three different formulas. Fortunately, one of them is correct!
    σ2 = σ^2
    c′ = f(c) + ddf(c) * σ2 / 2
    σ′² = df(c)^2 * σ2 + ddf(c)^2 * σ2^2 / 2
    return c′, sqrt(σ′²)
end

abs(a::GVar) = GVar(abs(a.center), a.σ)
abs2(a::GVar) = a^2

min(a::GVar, b::GVar) = lohi(GVar, min(loval(a), loval(b)), min(hival(a), hival(b)))
max(a::GVar, b::GVar) = lohi(GVar, max(loval(a), loval(b)), max(hival(a), hival(b)))
min(a::Real, b::GVar) = lohi(GVar, min(a, loval(b)), min(a, hival(b)))
min(a::GVar, b::Real) = min(b, a)
max(a::Real, b::GVar) = lohi(GVar, max(a, loval(b)), max(a, hival(b)))
max(a::GVar, b::Real) = max(b, a)

function Base.exp(a::GVar{<:AbstractFloat})
    isempty(a) && return a
    return GVar(meanvar(exp, exp, exp, a.center, a.σ)...)
end
Base.exp(a::GVar) = exp(float(a))

function Base.log(a::GVar{<:AbstractFloat})
    isempty(a) && return a
    return GVar(meanvar(log, inv, x->-1/x^2, a.center, a.σ)...)
end

function Base.sqrt(a::GVar{<:AbstractFloat})
    isempty(a) && return a
    return GVar(meanvar(sqrt, x -> 1/(2*sqrt(x)), x -> -1/(4 * sqrt(x^3)), a.center, a.σ)...)
end

end # module
