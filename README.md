# GaussianRandomVariables

This implements a numeric type, `GVar`, representing a Gaussian (normal) random variable. Some elementary mathematical functions of these variables are implemented, and these also return `GVar`s giving the approximate distribution of the output.

Demo:

```julia
julia> using GaussianRandomVariables

julia> x = 5 ± 1
5.0 ± 1.0

julia> 1/x
0.20800000000000002 ± 0.041569219381653054

```

Related package: [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) is a much more fully-developed and featureful package which also handles arithmetic with Gaussian random variables. However, it implements first-order (linear) error propagation, which leads to different rules of arithmetic: compare

```julia
julia> using GaussianRandomVariables

julia> (0 ± 1)^2
1.0 ± 1.4142135623730951
```

with

```julia
julia> using Measurements

julia> (0 ± 1)^2
0.0 ± 0.0
```

Measurements is recommended for most users, but GaussianRandomVariables can be recommended if second-order accuracy matters in your application.

GaussianRandomVariables is built on top of [ThickNumbers](https://github.com/timholy/ThickNumbers.jl), and the API for working with `GVar`s is described in detail there.
