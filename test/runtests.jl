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

ispositive(x) = x > 0

@testset "GaussianRandomVariables.jl" begin
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
    @test rad(1/x)^2 ≈ rad(x)^2 / mid(x)^4 + 2 * rad(x)^4 / mid(x)^6
    @test sqrt(x) ⩪ exp(0.5 * log(x)) rtol=1e-2
    @test sqrt(x) ⩪ x^0.5

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
