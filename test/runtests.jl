using GaussianRandomVariables
using ThickNumbers
using Test

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
end
