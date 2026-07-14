module GaussianRandomVariables

using ThickNumbers

import Base: +, -, *, /, //, ^, inv
import Base: abs, abs2, max, min, sqrt
import Base: log, exp

export GVar, ±
export skewness, moment_error, distrust

const BaseReals = Union{AbstractFloat, Integer, AbstractIrrational, Rational}

## Type and constructors

"""
    GVar(center, σ)
    center ± σ

A Gaussian random variable with mean `center` and standard deviation `σ`.

Arithmetic on `GVar`s propagates the distribution of the result to second order,
so `(0 ± 1)^2` is `1 ± √2` rather than the first-order answer `0 ± 0`.

A `GVar` additionally carries two *reliability diagnostics*, which accumulate
through arithmetic and say how far the result has drifted from being a genuine
Gaussian:

- [`skewness`](@ref), the standardized third cumulant. Zero for a Gaussian.
- [`moment_error`](@ref), an estimate of the absolute error in `center`.

[`distrust`](@ref) combines them into a single dimensionless number. When it is
large, `center` and `σ` do not describe the true distribution and should not be
believed.
"""
struct GVar{T<:Real} <: ThickNumber{T}
    center::T
    σ::T
    κ3::T      # third cumulant, equal to the third central moment
    err::T     # estimated error in `center`, accumulated from Taylor truncation

    function GVar{T}(center, σ, κ3, err) where T<:Real
        new{T}(center, σ, κ3, err)
    end
end

GVar{T}(center, σ) where T<:Real = GVar{T}(center, σ, zero(T), zero(T))
GVar{T}(x::Real) where T<:Real = GVar{T}(x, zero(T))

GVar(center::T, σ::T, κ3::T, err::T) where T<:Real = GVar{T}(center, σ, κ3, err)
GVar(center::Real, σ::Real, κ3::Real, err::Real) = GVar(promote(center, σ, κ3, err)...)
GVar(center::T, σ::T) where T<:Real = GVar{T}(center, σ)
GVar(center::Real, σ::Real) = GVar(promote(center, σ)...)
GVar(x::T) where T<:Real = GVar{T}(x)

# Integer and irrational parameters float-promote: the moment formulas divide.
GVar(center::T, σ::T) where T<:Integer = GVar(float(center), float(σ))
GVar(center::T, σ::T) where T<:Irrational = GVar(float(center), float(σ))
GVar(c::T, s::T, k::T, e::T) where T<:Integer = GVar(float(c), float(s), float(k), float(e))
GVar(c::T, s::T, k::T, e::T) where T<:Irrational = GVar(float(c), float(s), float(k), float(e))

GVar(x::GVar) = x
GVar{T}(x::GVar) where T = GVar{T}(x.center, x.σ, x.κ3, x.err)

±(center::Real, σ::Real) = GVar(center, σ)

Base.promote_rule(::Type{GVar{T}}, ::Type{GVar{S}}) where {T<:Real,S<:Real} = GVar{promote_type(T,S)}
Base.promote_rule(::Type{GVar{T}}, ::Type{S}) where {T<:Real,S<:BaseReals} = GVar{promote_type(T,S)}

Base.AbstractFloat(a::GVar{<:AbstractFloat}) = a
Base.AbstractFloat(a::GVar) = GVar(AbstractFloat(a.center), AbstractFloat(a.σ),
                                   AbstractFloat(a.κ3), AbstractFloat(a.err))

## Reliability diagnostics

"""
    skewness(a::GVar)

The standardized third cumulant of `a`, `κ3/σ^3`. It is zero for a true Gaussian
and grows as arithmetic drives the distribution away from one.

A nonzero skewness means the span `mid(a) ± k*rad(a)` is misplaced: to leading
order the `k`-sigma quantile sits at `mid(a) + rad(a)*(k + (k^2-1)*skewness(a)/6)`.
A large negative skewness is a long tail toward *low* values, i.e. the objective
may dip far below `loval(a)`.
"""
function skewness(a::GVar{T}) where T
    iszero(a.κ3) && return zero(float(T))
    return a.κ3 / a.σ^3
end

"""
    moment_error(a::GVar)

An estimate of the absolute error in `mid(a)`, accumulated across the operations
that produced `a`.

Arithmetic on `GVar`s is exact when each operation is at most quadratic over the
spread of its input, so this tracks the leading *neglected* Taylor term,
`|f''''(c)|σ^4/8`, propagated forward by `|f'(c)|` at each subsequent step.
"""
moment_error(a::GVar) = a.err

"""
    distrust(a::GVar)

A dimensionless measure of how much the Gaussian parameters of `a` can be
believed. It is zero when `a` is exactly Gaussian and grows without bound as the
description degrades; **smaller is better**.

Roughly, `distrust` is the displacement of `a`'s quantiles measured in units of
`rad(a)`: values ≲0.01 mean `mid(a)` and `rad(a)` are solid, while values ≳1 mean
they are fiction and the domain must be subdivided regardless of what they say.

It sums the two defects, which are also available on their own:
`moment_error(a)/rad(a)` (Taylor truncation) and `abs(skewness(a))/6` (the
leading Cornish–Fisher displacement from non-Gaussianity).

`distrust` is a heuristic, not a bound. It cannot see the error caused by
treating repeated appearances of a variable as independent, as in `x*x`.
"""
function distrust(a::GVar{T}) where T
    iszero(a.err) && iszero(a.κ3) && return zero(float(T))
    d = a.err / a.σ + abs(a.κ3) / (6 * a.σ^3)
    # A NaN here means the moments carry no information at all (e.g. `σ` is
    # infinite); that is the least trustworthy state, not an unremarkable one.
    return isnan(d) ? convert(float(T), Inf) : d
end

## ThickNumbers API

ThickNumbers.loval(a::GVar) = a.center - a.σ
ThickNumbers.hival(a::GVar) = a.center + a.σ
ThickNumbers.lohi(::Type{G}, lo, hi) where G<:GVar = G((lo + hi)/2, nextfloat((hi - lo)/2))
ThickNumbers.midrad(::Type{G}, center, σ) where G<:GVar = G(center, σ)
ThickNumbers.basetype(::Type{GVar{T}}) where T = GVar
ThickNumbers.basetype(::Type{GVar}) = GVar
ThickNumbers.emptyset(::Type{G}) where G<:GVar = G(0, -1)

# ThickNumbers' generic `hull`, `intersect`, `typemin` and `typemax` construct a
# result with the two-argument call `TN(lo, hi)`, but `GVar`'s two-argument
# constructor is `GVar(center, σ)`. Build these from the span directly.
Base.typemin(::Type{GVar{T}}) where T<:Real = GVar{T}(typemin(T), zero(T))
Base.typemax(::Type{GVar{T}}) where T<:Real = GVar{T}(typemax(T), zero(T))

# A hull or intersection is a set operation, not a Gaussian: it has no meaningful
# skew, and it must not gain reliability that neither operand had.
_span(::Type{G}, lo, hi, err) where G<:GVar = G((lo + hi)/2, (hi - lo)/2, zero(lo), err)

function Base.intersect(a::GVar{T}, b::GVar{T}) where T
    isdisjoint(a, b) && return emptyset(GVar{T})
    return _span(GVar{T}, max(loval(a), loval(b)), min(hival(a), hival(b)), max(a.err, b.err))
end
ThickNumbers.hull(a::GVar{T}, b::GVar{T}) where T =
    _span(GVar{T}, min(loval(a), loval(b)), max(hival(a), hival(b)), max(a.err, b.err))

## Display

function Base.show(io::IO, a::GVar)
    print(io, a.center, " ± ", a.σ)
    d = distrust(a)
    iszero(d) || print(io, " (distrust ", round(d; sigdigits=2), ")")
end

## Trait functions and constants

Base.zero(::GVar{T}) where T<:Real = zero(GVar{T})
Base.zero(::Type{GVar{T}}) where T<:Real = GVar(zero(T), zero(T))
Base.oneunit(::GVar{T}) where T<:Real = oneunit(GVar{T})
Base.oneunit(::Type{GVar{T}}) where T<:Real = GVar(oneunit(T), zero(T))

Base.real(a::GVar) = a
Base.conj(a::GVar) = a

function Base.hash(x::GVar, h::UInt)
    magic = Int === Int64 ? 0x21f1b1afbb07c31a : 0x237f8305
    h += magic
    return hash(x.center, hash(x.σ, hash(x.κ3, hash(x.err, h))))
end

## Moment propagation

# For x ~ N(c, σ²) with third cumulant κ3, and smooth f, ordering terms by h where
# σ ~ h and κ3 ~ h³:
#
#   E[f]   = f + f''σ²/2 + f'''κ3/6                     + O(h⁴)
#   Var[f] = f'²σ² + f'f''κ3 + f''²σ⁴/2 + f'f'''σ⁴      + O(h⁵)
#   κ3[f]  = f'³κ3 + 3f'²f''σ⁴ + f''³σ⁶                 + O(h⁷)
#
# The leading neglected term of the mean, f''''σ⁴/8, is not added to the result;
# it is accumulated into `err` as the error estimate. Both σ⁴ terms of the
# variance are the same order, so keeping one without the other understates σ.
function gmap(f::F, df::D1, ddf::D2, dddf::D3, d4f::D4, a::GVar) where {F,D1,D2,D3,D4}
    c, σ, κ3, err = a.center, a.σ, a.κ3, a.err
    f1, f2, f3, f4 = df(c), ddf(c), dddf(c), d4f(c)
    σ2 = σ^2
    return assemble(f(c) + f2*σ2/2 + f3*κ3/6,
                    f1^2*σ2 + f1*f2*κ3 + f2^2*σ2^2/2 + f1*f3*σ2^2,
                    f1^3*κ3 + 3*f1^2*f2*σ2^2 + f2^3*σ2^3,
                    abs(f1)*err + abs(f4)*σ2^2/8)
end

# A negative variance means the expansion has broken down completely, which
# happens only when the input is far from Gaussian. Surrender the moments rather
# than report a plausible-looking answer.
function assemble(m, v, κ3, err)
    v < zero(v) && return GVar(m, zero(v), zero(κ3), oftype(err, Inf))
    return GVar(m, sqrt(v), κ3, err)
end

# E[(c + δ)^n] for δ ~ N(0, σ²). Odd powers of δ vanish and E[δ^k] = σ^k (k-1)!!,
# so the sum terminates: for polynomial f the Gaussian moments are exact.
function gaussrawmoment(n::Integer, c::T, σ::T) where T<:AbstractFloat
    s = zero(T)
    binom = one(T)         # C(n, k)
    dfact = one(T)         # (k-1)!!
    for k in 0:2:n
        if k > 0
            binom *= T(n - k + 2) * T(n - k + 1) / (T(k) * T(k - 1))
            dfact *= T(k - 1)
        end
        s += binom * c^(n - k) * σ^k * dfact
    end
    return s
end

## Addition and subtraction

+(a::GVar{<:Real}, b::Real) = GVar(a.center + b, a.σ, a.κ3, a.err)
+(a::GVar{<:Integer}, b::Integer) = float(a) + float(b)
+(b::Real, a::GVar{<:Real}) = a + b

-(a::GVar{<:Real}) = GVar(-a.center, a.σ, -a.κ3, a.err)
-(a::GVar{<:Real}, b::Real) = GVar(a.center - b, a.σ, a.κ3, a.err)
-(a::GVar{<:Integer}, b::Integer) = float(a) - float(b)
-(b::Real, a::GVar{<:Real}) = GVar(b - a.center, a.σ, -a.κ3, a.err)
-(b::Integer, a::GVar{<:Integer}) = float(b) - float(a)

# Cumulants of independent variables add.
function +(a::GVar{<:Real}, b::GVar{<:Real})
    @fastmath begin
        ret = GVar(a.center + b.center, sqrt(a.σ^2 + b.σ^2), a.κ3 + b.κ3, a.err + b.err)
        return (isempty(a) | isempty(b)) ? emptyset(ret) : ret
    end
end
+(a::GVar{<:Integer}, b::GVar{<:Integer}) = float(a) + float(b)

function -(a::GVar{<:Real}, b::GVar{<:Real})
    @fastmath begin
        ret = GVar(a.center - b.center, sqrt(a.σ^2 + b.σ^2), a.κ3 - b.κ3, a.err + b.err)
        return (isempty(a) | isempty(b)) ? emptyset(ret) : ret
    end
end
-(a::GVar{<:Integer}, b::GVar{<:Integer}) = float(a) - float(b)

## Multiplication

# Restrict to BaseReals so we avoid ambiguities with ForwardDiff
*(x::BaseReals, a::GVar{<:Real}) = GVar(x*a.center, abs(x)*a.σ, x^3*a.κ3, abs(x)*a.err)
*(a::GVar{<:Real}, x::BaseReals) = x*a
# Prevent overflow by promoting to float
*(x::Integer, a::GVar{<:Integer}) = float(x)*float(a)
*(a::GVar{<:Integer}, x::Integer) = float(x)*float(a)

# For independent a and b these moments are exact given the first three cumulants
# of each: no fourth moment appears.
function *(a::GVar{<:Real}, b::GVar{<:Real})
    @fastmath begin
        ca, σa, κa, ea = a.center, a.σ, a.κ3, a.err
        cb, σb, κb, eb = b.center, b.σ, b.κ3, b.err
        v = ca^2*σb^2 + cb^2*σa^2 + σa^2*σb^2
        κ3 = ca^3*κb + cb^3*κa + κa*κb + 3*ca*σa^2*κb + 3*cb*σb^2*κa + 6*ca*cb*σa^2*σb^2
        ret = GVar(ca*cb, sqrt(v), κ3, abs(cb)*ea + abs(ca)*eb + ea*eb)
        return (isempty(a) | isempty(b)) ? emptyset(ret) : ret
    end
end
*(a::GVar{<:Integer}, b::GVar{<:Integer}) = float(a)*float(b)   # prevent overflow

## Division

/(a::GVar{<:Real}, x::Real) = GVar(a.center/x, a.σ/abs(x), a.κ3/x^3, a.err/abs(x))

function inv(a::GVar{T}) where T<:AbstractFloat  # must be AbstractFloat so typemax gives Inf
    z = zero(1/oneunit(T))
    isempty(a) && return emptyset(basetype(a){typeof(z)})
    # Straddling zero, the reciprocal has no finite moments at all.
    zero(T) ∈ a && return GVar(iszero(a.center) ? z : 1/a.center, typemax(T), zero(T), typemax(T))
    return gmap(inv, x -> -1/x^2, x -> 2/x^3, x -> -6/x^4, x -> 24/x^5, a)
end
inv(a::GVar{<:Real}) = inv(float(a))

/(a::Real, b::GVar{<:Real}) = a*inv(b)
/(a::GVar{<:Real}, b::GVar{<:Real}) = a*inv(b)
//(a::GVar, b::GVar) = a / b    # to deal with rationals

## Powers

Base.literal_pow(::typeof(^), x::GVar, ::Val{0}) = oneunit(x)
Base.literal_pow(::typeof(^), x::GVar, ::Val{1}) = x
Base.literal_pow(::typeof(^), x::GVar, ::Val{p}) where p = x^p

function ^(a::GVar{T}, p::Integer) where T<:AbstractFloat
    isempty(a) && return a
    p == 0 && return oneunit(a)
    p == 1 && return a
    p < 0 && return inv(a^(-p))
    3p <= 60 || throw(DomainError(p, "exponents above 20 overflow the moment sum; rescale or use a float exponent"))
    c, σ, κ3, err = a.center, a.σ, a.κ3, a.err
    # x^p is a polynomial, so the Gaussian moments terminate and are exact. The
    # only error is the input's own, plus the leading correction for its skew.
    m = gaussrawmoment(p, c, σ)
    m2 = gaussrawmoment(2p, c, σ)
    m3 = gaussrawmoment(3p, c, σ)
    f1 = p*c^(p-1)
    f2 = p*(p-1)*c^(p-2)
    f3 = p >= 3 ? p*(p-1)*(p-2)*c^(p-3) : zero(T)
    return assemble(m + f3*κ3/6,
                    (m2 - m^2) + f1*f2*κ3,
                    (m3 - 3*m*m2 + 2*m^3) + f1^3*κ3,
                    abs(f1)*err)
end
^(a::GVar{<:Real}, p::Integer) = float(a)^p

function ^(a::GVar, p::Real)
    isinteger(p) && return a^Int(p)
    isempty(a) && return a
    return gmap(x -> x^p,
                x -> p*x^(p-1),
                x -> p*(p-1)*x^(p-2),
                x -> p*(p-1)*(p-2)*x^(p-3),
                x -> p*(p-1)*(p-2)*(p-3)*x^(p-4), a)
end

## Functions

function abs(a::GVar)
    isempty(a) && return a
    # Away from zero `abs` merely reflects. Straddling zero the result is a folded
    # normal, which is not Gaussian at all; the center is then wrong by O(σ).
    straddles = zero(a.center) ∈ a
    err = straddles ? a.err + a.σ : a.err
    return GVar(abs(a.center), a.σ, sign(a.center)*a.κ3, err)
end
abs2(a::GVar) = a^2

# `min` and `max` of Gaussians are not Gaussian. Keep the span, report no skew,
# and inherit the larger center error.
for (f, lo, hi) in ((:min, :min, :min), (:max, :max, :max))
    @eval begin
        function $f(a::GVar, b::GVar)
            r = lohi(GVar, $lo(loval(a), loval(b)), $hi(hival(a), hival(b)))
            return GVar(r.center, r.σ, zero(r.κ3), max(a.err, b.err))
        end
        function $f(a::Real, b::GVar)
            r = lohi(GVar, $lo(a, loval(b)), $hi(a, hival(b)))
            return GVar(r.center, r.σ, zero(r.κ3), b.err)
        end
        $f(a::GVar, b::Real) = $f(b, a)
    end
end

function exp(a::GVar{<:AbstractFloat})
    isempty(a) && return a
    c, σ, κ3, err = a.center, a.σ, a.κ3, a.err
    σ2 = σ^2
    # The exponential of a Gaussian is exactly lognormal, so there is no
    # truncation in σ: only the input's skew and center error propagate.
    e2 = exp(σ2)
    m = exp(c + σ2/2)
    v = (e2 - 1) * exp(2c + σ2)
    k = (e2 - 1)^2 * (e2 + 2) * exp(3c + 3σ2/2)
    f1 = exp(c)     # f' = f'' = f''' for exp
    return assemble(m + f1*κ3/6, v + f1^2*κ3, k + f1^3*κ3, m*err)
end

log(a::GVar{<:AbstractFloat}) =
    isempty(a) ? a : gmap(log, x -> 1/x, x -> -1/x^2, x -> 2/x^3, x -> -6/x^4, a)

sqrt(a::GVar{<:AbstractFloat}) =
    isempty(a) ? a : gmap(sqrt, x -> 1/(2*sqrt(x)), x -> -1/(4*sqrt(x^3)),
                          x -> 3/(8*sqrt(x^5)), x -> -15/(16*sqrt(x^7)), a)

end # module
