using GaussianRandomVariables
using ThickNumbers
using Statistics
using HypothesisTests
using StableRNGs
using Test

# These comparisons are Monte Carlo, so an unseeded generator turns the suite into a
# coin flip. Julia's own streams are not reproducible across releases, so seed a
# StableRNG to keep a given commit's result identical on every version.
const rng = StableRNG(20240614)

# The loop below runs roughly 70 hypothesis tests, so a per-test threshold of 0.001
# rejects a correct implementation in about 7% of runs. 1e-5 is the corresponding
# Bonferroni level; it still detects a 25% error in σ on every draw.
function testscalar(f, μ, σ; n=1000, pval = 1e-5, filter=Returns(true), rtol=nothing)
    x = Base.filter(filter, μ .+ σ .* randn(rng, n))
    @assert length(x) >= 0.9 * n
    y = f.(x)
    g = GVar(μ, σ)
    fg = f(g)
    if rtol === nothing
        z = mid(fg) .+ rad(fg) .* randn(rng, n)
        return pvalue(EqualVarianceTTest(y, z)) > pval && pvalue(LeveneTest(y, z)) > pval
    end
    return isapprox(mean(y), mid(fg); rtol) && isapprox(std(y), rad(fg); rtol)
end

# Third central moment of a large sample, for checking `κ3` propagation.
function testskew(f, μ, σ; n=10^7, rtol=0.05, filter=Returns(true))
    x = Base.filter(filter, μ .+ σ .* randn(rng, n))
    y = f.(x)
    m = mean(y)
    return isapprox(mean((y .- m).^3), f(GVar(μ, σ)).κ3; rtol)
end

ispositive(x) = x > 0

@testset "GaussianRandomVariables.jl" begin
    @testset "arithmetic" begin
        x = 3 ± 1
        e = GVar(0, -1)          # empty
        @test isempty(e)
        @test isempty(e^1)
        @test isempty(e^2)
        @test isempty(^(e, 2))   # non-literal evaluation
        @test x^2 ⩪ ^(x, 2)      # non-literal
        @test mid(x*x) < mid(x^2)
        @test rad(x*x) < rad(x^2)
        @test rad(GVar(0, 1)^2) ≈ sqrt(2)
        @test mid(x - x) == 0
        @test rad(x - x) ≈ sqrt(2)
        @test mid(1/x) ≈ 1/mid(x) + rad(x)^2/mid(x)^3
        # Var[1/x] = σ²/c⁴ + 2σ⁴/c⁶ + f'f'''σ⁴, and f'f''' = 6/c⁶ for f = inv
        @test rad(1/x)^2 ≈ rad(x)^2/mid(x)^4 + 8*rad(x)^4/mid(x)^6
        @test sqrt(x) ⩪ exp(0.5 * log(x)) rtol=1e-2
        @test sqrt(x) ⩪ x^0.5
        @test mid(-x) == -mid(x) && rad(-x) == rad(x)
        @test x^0 ⩪ oneunit(x)
        @test x^1 ⩪ x
        @test 1/(1/x) ⩪ x rtol=0.1
    end

    @testset "distributions" begin
        for μ in (2, 3, 10), σ in (0.1, 0.25, 1)
            # `GVar` propagates moments to second order, which holds only while the spread
            # is small compared to the distance from `μ` to a singularity or a sign change.
            # At μ = 3σ the sampled distribution of `log(x)` already has a tail reaching
            # tens of radii below `mid`, and these tests rightly reject it.
            μ >= 4σ || continue
            @test testscalar(x -> x^2, μ, σ)
            @test testscalar(x -> x^1.8, μ, σ; filter=ispositive)
            @test testscalar(exp, μ, σ; rtol=0.03, n=10^6)
            @test testscalar(log, μ, σ; filter=ispositive)
            @test testscalar(sqrt, μ, σ; filter=ispositive)
        end
    end

    # x^p is a polynomial, so the Gaussian moments terminate: mean, variance and
    # third cumulant are all exact.
    @testset "exact moments of integer powers" begin
        c, σ = 3.0, 0.75
        y = GVar(c, σ)^2
        @test mid(y) ≈ c^2 + σ^2
        @test rad(y)^2 ≈ 4c^2*σ^2 + 2σ^4
        @test y.κ3 ≈ 24c^2*σ^4 + 8σ^6
        @test moment_error(y) == 0                 # nothing is neglected
        # x² of a zero-mean Gaussian is σ²·χ²₁, whose skewness is √8
        @test skewness(GVar(0, 1)^2) ≈ sqrt(8)
        @test testskew(x -> x^3, 2.0, 0.5)
    end

    @testset "reliability diagnostics" begin
        # A box straight from a bisection is a plain Gaussian: nothing to distrust.
        b = lohi(GVar, 0.0, 1.0)
        @test distrust(b) == 0
        @test skewness(b) == 0
        @test moment_error(b) == 0

        # Quadratics are handled exactly, so they must not be flagged...
        @test moment_error((3 ± 0.5)^2 + 1) == 0
        # ...while a strongly-curved function over a wide box must be.
        @test distrust(exp(1 ± 1.5)) > 1
        @test distrust(sqrt(3 ± 0.01)) < 0.01

        # distrust decreases monotonically as the domain is subdivided, which is
        # what makes "bisect until trustworthy" terminate.
        ds = [distrust(1 / ((lohi(GVar, 0.5 - w, 0.5 + w) - 0.2)^2 + 0.05)) for w in (0.4, 0.2, 0.1, 0.05, 0.025)]
        @test issorted(ds; rev=true)
        @test last(ds) < 0.1

        # exp of a Gaussian is exactly lognormal: the moments are right, but the
        # result is badly skewed, so mid ∓ rad is not a sensible span.
        y = exp(1 ± 0.7)
        @test moment_error(y) == 0
        @test skewness(y) > 2
        @test distrust(y) > 0.4

        # Straddling zero, 1/x has no finite moments at all.
        @test distrust(1 / (0 ± 1)) == Inf
    end

    @testset "set operations" begin
        a, b = GVar(3.0, 1.0), GVar(3.5, 1.0)   # spans [2,4] and [2.5,4.5]
        @test loval(hull(a, b)) ≈ 2.0 && hival(hull(a, b)) ≈ 4.5
        @test loval(intersect(a, b)) ≈ 2.5 && hival(intersect(a, b)) ≈ 4.0
        @test isempty(intersect(GVar(0.0, 1.0), GVar(10.0, 1.0)))
        @test loval(typemax(GVar{Float64})) == Inf
        @test hival(typemin(GVar{Float64})) == -Inf
        # A result built from the span alone would report the diagnostics of a
        # fresh box; combining values must not manufacture reliability.
        u = 1 / ((3 ± 0.9)^2 - 8)     # moment_error > 0
        @test moment_error(u) > 0
        @test moment_error(hull(u, GVar(1.0, 1.0))) >= moment_error(u)
        @test moment_error(intersect(u, GVar(mid(u), 10rad(u)))) >= moment_error(u)
    end

    # Cumulants of independent variables add.
    @testset "sums and differences" begin
        a, b = GVar(1.0, 0.5), GVar(-2.0, 1.2)
        @test mid(a + b) ≈ -1.0
        @test rad(a + b) ≈ sqrt(0.5^2 + 1.2^2)
        @test mid(a - b) ≈ 3.0
        @test rad(a - b) ≈ sqrt(0.5^2 + 1.2^2)
        @test @inferred(a + b) isa GVar{Float64}

        # Skew adds under `+` and subtracts under `-`; so does the center error.
        s = (2 ± 0.5)^3          # κ3 > 0, err == 0 (a polynomial is exact)
        t = exp(1 ± 0.3)         # κ3 > 0
        @test (s + t).κ3 ≈ s.κ3 + t.κ3
        @test (s - t).κ3 ≈ s.κ3 - t.κ3
        u = sqrt(4 ± 0.5)        # err > 0
        @test moment_error(u) > 0
        @test moment_error(u + s) ≈ moment_error(u) + moment_error(s)

        # Sampled independently, the sum of two Gaussians is Gaussian.
        x = 1.0 .+ 0.5 .* randn(rng, 1000)
        y = -2.0 .+ 1.2 .* randn(rng, 1000)
        z = mid(a + b) .+ rad(a + b) .* randn(rng, 1000)
        @test pvalue(EqualVarianceTTest(x .+ y, z)) > 1e-5
        @test pvalue(LeveneTest(x .+ y, z)) > 1e-5

        e = GVar(0.0, -1.0)
        @test isempty(e + a) && isempty(a + e)
        @test isempty(e - a) && isempty(a - e)
    end

    @testset "abs" begin
        # Away from zero, `abs` merely reflects.
        @test abs(GVar(-3.0, 0.1)) ⩪ GVar(3.0, 0.1)
        @test abs((-2 ± 0.05)^3).κ3 ≈ -((-2 ± 0.05)^3).κ3   # reflection flips the skew

        # Straddling zero the result is a folded normal, whose moments are exact.
        # For a zero-mean input this is the half-normal: mean σ√(2/π), variance
        # σ²(1 - 2/π), skewness √2(4 - π)/(π - 2)^(3/2).
        h = abs(0 ± 1)
        @test mid(h) ≈ sqrt(2/π)
        @test rad(h) ≈ sqrt(1 - 2/π)
        @test skewness(h) ≈ sqrt(2)*(4 - π)/(π - 2)^(3//2)
        # The center is now right, so nothing is charged to `err`; only the shape
        # is non-Gaussian, and `distrust` sees that through the skew alone.
        @test moment_error(h) == 0
        @test distrust(h) ≈ abs(skewness(h))/6

        @test testscalar(abs, 0.0, 1.0; rtol=0.02, n=10^6)
        @test testscalar(abs, 0.5, 1.0; rtol=0.02, n=10^6)
        @test testskew(abs, -0.7, 1.0)

        # Reflection and folding must agree where they meet.
        @test abs(GVar(8.0, 1.0)) ⩪ abs(GVar(prevfloat(8.0), 1.0))
        @test abs(GVar(0.0, 0.0)) ⩪ GVar(0.0, 0.0)      # σ == 0: no fold to compute
        @test isempty(abs(GVar(0.0, -1.0)))
        @test abs2(3 ± 0.5) ⩪ (3 ± 0.5)^2
    end

    # sin and cos of a Gaussian are not Gaussian, but their first three moments are
    # exact, so mid, rad and κ3 match closed forms and nothing is charged to `err`.
    @testset "sin and cos" begin
        c, σ = 0.7, 0.5
        d1, d2 = exp(-σ^2/2), exp(-2σ^2)
        s = sin(GVar(c, σ))
        @test mid(s) ≈ sin(c)*d1
        @test rad(s)^2 ≈ (1 - cos(2c)*d2)/2 - (sin(c)*d1)^2
        @test moment_error(s) == 0
        k = cos(GVar(c, σ))
        @test mid(k) ≈ cos(c)*d1
        @test rad(k)^2 ≈ (1 + cos(2c)*d2)/2 - (cos(c)*d1)^2
        @test moment_error(k) == 0

        # A zero-mean input is symmetric: sin has zero mean and no skew, while cos
        # peaks at 1 and folds its mass downward, so it is skewed toward low values.
        @test mid(sin(0 ± σ)) == 0
        @test sin(0 ± σ).κ3 == 0
        @test mid(cos(0 ± σ)) ≈ exp(-σ^2/2)
        @test skewness(cos(0 ± σ)) < 0

        # A vanishing spread is an ordinary number.
        @test sin(GVar(c, 0.0)) ⩪ GVar(sin(c), 0.0)
        @test cos(GVar(c, 0.0)) ⩪ GVar(cos(c), 0.0)

        # `sincos` returns the pair; the empty set propagates through all three.
        sc = sincos(GVar(c, σ))
        @test sc[1] ⩪ sin(GVar(c, σ)) && sc[2] ⩪ cos(GVar(c, σ))
        e = GVar(0.0, -1.0)
        @test isempty(sin(e)) && isempty(cos(e))
        @test all(isempty, sincos(e))

        @test @inferred(sin(GVar(1.0, 0.5))) isa GVar{Float64}
        @test @inferred(sincos(GVar(1.0, 0.5))) isa Tuple{GVar{Float64},GVar{Float64}}

        # Sampled moments confirm the closed forms across a range of spreads.
        for (μ, σ) in ((0.6, 0.2), (1.0, 0.5), (2.0, 0.8), (0.9, 1.0))
            @test testscalar(sin, μ, σ; rtol=0.02, n=10^6)
            @test testscalar(cos, μ, σ; rtol=0.02, n=10^6)
        end
        for (μ, σ) in ((1.0, 0.5), (2.0, 0.8), (0.9, 1.0))
            @test testskew(sin, μ, σ)
            @test testskew(cos, μ, σ)
        end
    end

    @testset "min and max" begin
        a, b = GVar(3.0, 1.0), GVar(3.5, 1.0)    # spans [2,4] and [2.5,4.5]
        @test loval(min(a, b)) ≈ 2.0 && hival(min(a, b)) ≈ 4.0
        @test loval(max(a, b)) ≈ 2.5 && hival(max(a, b)) ≈ 4.5
        # A min/max is a set operation, not a Gaussian: no skew, and the larger
        # center error survives.
        u = 1 / ((3 ± 0.9)^2 - 8)
        @test moment_error(u) > 0
        @test min(u, a).κ3 == max(u, a).κ3 == 0
        @test moment_error(min(u, a)) == moment_error(max(u, a)) == moment_error(u)

        # Mixed with a plain number, in either argument order.
        @test min(2.0, a) ⩪ min(a, 2.0)
        @test max(2.0, a) ⩪ max(a, 2.0)
        @test loval(max(2.0, a)) ≈ 2.0 && hival(max(2.0, a)) ≈ 4.0
        # A degenerate span has zero radius: `min(a, 2.0)` is exactly 2.0, since
        # every value in `a` is at least 2.
        @test min(a, 2.0) ⩪ GVar(2.0, 0.0)
        @test rad(min(a, 2.0)) == 0
    end

    # `lohi` must return a span that contains [lo, hi] even after rounding, without
    # inflating one that is already exact.
    @testset "lohi rounding" begin
        @test rad(lohi(GVar, 1.0, 1.0)) == 0
        @test lohi(GVar, 0.0, 1.0) ⩪ GVar(0.5, 0.5)
        @test rad(lohi(GVar, 0.0, 1.0)) == 0.5
        for (lo, hi) in ((0.0, 1.0), (-1.0, 1e300), (1e100, nextfloat(1e100)),
                         (1e-300, 1.0), (-2.0, -1.0), (0.1, 0.1 + 1e-17))
            g = lohi(GVar, lo, hi)
            @test loval(g) <= lo && hival(g) >= hi
        end
    end

    @testset "scalar arithmetic and promotion" begin
        x = 3.0 ± 1.0
        @test x + 2 ⩪ 2 + x
        @test mid(x + 2) == 5.0 && rad(x + 2) == 1.0
        @test x - 2 ⩪ GVar(1.0, 1.0)
        @test 2 - x ⩪ GVar(-1.0, 1.0)
        @test (2 - x).κ3 == -x.κ3
        @test x * 2 ⩪ 2 * x
        @test mid(x * 2) == 6.0 && rad(x * 2) == 2.0
        @test x / 2 ⩪ 0.5 * x
        @test (-2 * (2 ± 0.5)^3).κ3 ≈ -8 * ((2 ± 0.5)^3).κ3   # κ3 scales as the cube
        @test x / (2 ± 0.5) ⩪ x * inv(2 ± 0.5)
        @test x // (2 ± 0.5) ⩪ x / (2 ± 0.5)

        # `GVar(1, 2)` already float-promotes, so an integer-backed `GVar` can only
        # come from the parametric constructor. Arithmetic on one promotes rather
        # than overflowing: the moment formulas divide, and κ3 grows as the cube.
        i, j = GVar{Int}(1, 2), GVar{Int}(3, 4)
        @test i isa GVar{Int}
        @test @inferred(i + j) ⩪ GVar(4.0, sqrt(20.0))
        @test @inferred(i - j) ⩪ GVar(-2.0, sqrt(20.0))
        @test @inferred(i * j) isa GVar{Float64}
        @test i + 2 ⩪ GVar(3.0, 2.0)
        @test i - 2 ⩪ GVar(-1.0, 2.0)
        @test 2 - i ⩪ GVar(1.0, 2.0)
        @test 2 * i ⩪ i * 2 ⩪ GVar(2.0, 4.0)
        @test i * 2 isa GVar{Float64}
        @test inv(GVar{Int}(2, 1)) ⩪ inv(GVar(2.0, 1.0))
        @test GVar{Int}(3, 1)^2 ⩪ GVar(3.0, 1.0)^2
        @test abs(GVar{Int}(-3, 1)) ⩪ abs(GVar(-3.0, 1.0))
        @test_throws "exponents above 20 overflow" GVar(1.0, 0.1)^21
    end

    @testset "display" begin
        # A trustworthy value shows only its span; a distrusted one says so.
        @test repr(GVar(1.0, 2.0)) == "1.0 ± 2.0"
        @test repr(GVar(1.0, 2.0, 4.8, 0.0)) == "1.0 ± 2.0 (distrust 0.1)"   # |κ3|/6σ³
        @test occursin("distrust", repr(exp(1 ± 1.5)))
    end

    @testset "traits" begin
        @test zero(GVar(1.0, 2.0)) === zero(GVar{Float64}) === GVar(0.0, 0.0)
        @test oneunit(GVar(1.0, 2.0)) === oneunit(GVar{Float64}) === GVar(1.0, 0.0)
        @test real(GVar(1.0, 2.0)) === GVar(1.0, 2.0)
        @test conj(GVar(1.0, 2.0)) === GVar(1.0, 2.0)
        @test basetype(GVar{Float64}) === basetype(GVar) === GVar
        @test promote_type(GVar{Float64}, GVar{Float32}) === GVar{Float64}
        @test promote_type(GVar{Float64}, Int) === GVar{Float64}
        @test AbstractFloat(GVar(1.0, 2.0)) === GVar(1.0, 2.0)
        @test AbstractFloat(GVar{Int}(1, 2)) === GVar(1.0, 2.0)

        # `hash` must see every field: two `GVar`s with the same span can still
        # differ in their diagnostics.
        @test hash(GVar(1.0, 2.0)) == hash(GVar(1.0, 2.0))
        @test hash(GVar(1.0, 2.0, 3.0, 0.0)) != hash(GVar(1.0, 2.0, 0.0, 0.0))
        @test hash(GVar(1.0, 2.0, 0.0, 3.0)) != hash(GVar(1.0, 2.0, 0.0, 0.0))
    end

    @testset "constructors" begin
        @test GVar(1, 2.0, 3, 4) === GVar(1.0, 2.0, 3.0, 4.0)
        @test GVar(1, 2, 3, 4) === GVar(1.0, 2.0, 3.0, 4.0)
        @test GVar(π, π, π, π) isa GVar{Float64}
        @test GVar(3.0) === GVar(3.0, 0.0)
        @test GVar(GVar(1.0, 2.0)) === GVar(1.0, 2.0)
        @test GVar{Float64}(1, 0) === GVar(1.0, 0.0)
        @test GVar(1, 2) === GVar(1.0, 2.0)
        @test GVar(π, π) isa GVar{Float64}
        @test GVar{Float32}(GVar(1.0, 2.0)) isa GVar{Float32}
        @test valuetype(GVar(1.0, 2.0)) === Float64
        @test convert(GVar{Float64}, 3) ⩪ GVar(3.0, 0.0)
    end
end
