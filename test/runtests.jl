using GaussianRandomVariables
using ThickNumbers
using HypothesisTests
using Test

function testscalar(f, μ, σ; n=1000, pval = 0.001)
    x = μ .+ σ .* randn(n)
    y = f.(x)
    g = GVar(μ, σ)
    fg = f(g)
    z = mid(fg) .+ rad(fg) .* randn(n)
    return pvalue(EqualVarianceTTest(y, z)) > pval && pvalue(LeveneTest(y, z)) > pval
end

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
    @test sqrt(x) ⩪ exp(0.5 * log(x)) rtol=1e-3
    @test sqrt(x) ⩪ x^0.5

    for μ in (2, 3, 10), σ in (0.25,)
        @test testscalar(x -> x^2, μ, σ)
        @test testscalar(x -> x^1.8, μ, σ)
        @test testscalar(exp, μ, σ)
        @test testscalar(log, μ, σ)
        @test testscalar(sqrt, μ, σ)
    end
end
