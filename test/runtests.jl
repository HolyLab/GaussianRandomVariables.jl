using GaussianRandomVariables
using ThickNumbers
using Statistics
using HypothesisTests
using Test

function testscalar(f, μ, σ; n=1000, pval = 0.001, filter=Returns(true), rtol=nothing)
    x = Base.filter(filter, μ .+ σ .* randn(n))
    @assert length(x) >= 0.9 * n
    y = f.(x)
    g = GVar(μ, σ)
    fg = f(g)
    if rtol === nothing
        z = mid(fg) .+ rad(fg) .* randn(n)
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
        @test testscalar(x -> x^2, μ, σ)
        @test testscalar(x -> x^1.8, μ, σ; filter=ispositive)
        @test testscalar(exp, μ, σ; rtol=0.03, n=10^6)
        μ - 3*σ >= 0 &&  @test testscalar(log, μ, σ; filter=ispositive)
        @test testscalar(sqrt, μ, σ; filter=ispositive)
    end
end
