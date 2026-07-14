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

    # ThickNumbers' generic `hull`/`intersect`/`typemin`/`typemax` build a result
    # with `TN(lo, hi)`, but `GVar(center, σ)` is a different parametrization.
    @testset "set operations" begin
        a, b = GVar(3.0, 1.0), GVar(3.5, 1.0)   # spans [2,4] and [2.5,4.5]
        @test loval(hull(a, b)) ≈ 2.0 && hival(hull(a, b)) ≈ 4.5
        @test loval(intersect(a, b)) ≈ 2.5 && hival(intersect(a, b)) ≈ 4.0
        @test isempty(intersect(GVar(0.0, 1.0), GVar(10.0, 1.0)))
        @test loval(typemax(GVar{Float64})) == Inf
        @test hival(typemin(GVar{Float64})) == -Inf
        # Combining two values must never manufacture reliability neither had.
        u = exp(1 ± 1.5)
        @test moment_error(hull(u, GVar(1.0, 1.0))) >= moment_error(u)
    end

    @testset "constructors" begin
        @test GVar{Float64}(1, 0) === GVar(1.0, 0.0)
        @test GVar(1, 2) === GVar(1.0, 2.0)
        @test GVar(π, π) isa GVar{Float64}
        @test GVar{Float32}(GVar(1.0, 2.0)) isa GVar{Float32}
        @test valuetype(GVar(1.0, 2.0)) === Float64
        @test convert(GVar{Float64}, 3) ⩪ GVar(3.0, 0.0)
    end
end
