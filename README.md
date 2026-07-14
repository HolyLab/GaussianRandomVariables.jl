# GaussianRandomVariables

This implements a numeric type, `GVar`, representing a Gaussian (normal) random variable. Some elementary mathematical functions of these variables are implemented, and these also return `GVar`s giving the approximate distribution of the output.

Demo:

```julia
julia> using GaussianRandomVariables

julia> x = 5 ± 1
5.0 ± 1.0

julia> 1/x
0.20800000000000002 ± 0.04595650117230423 (distrust 0.16)

```

## Reliability

Arithmetic on `GVar`s is exact whenever each operation is at most quadratic over the spread of its input, and whenever the input really is Gaussian. Neither holds in general, and when they fail the reported `center` and `σ` can be badly wrong, most dangerously by making the value look more tightly determined than it is.

A `GVar` therefore carries two diagnostics that accumulate through arithmetic:

- `skewness(a)`, the standardized third cumulant, is zero for a true Gaussian. It measures how far the result has drifted from the shape a `GVar` can represent. Its sign is informative: a large negative skewness is a long tail toward *low* values.
- `moment_error(a)` estimates the absolute error in `center`, by tracking the leading Taylor term each operation neglects.

`distrust(a)` combines them into one dimensionless number, roughly "how far `a`'s quantiles are displaced, in units of `rad(a)`". Smaller is better; zero means exactly Gaussian.

```julia
julia> distrust((3 ± 0.5)^2 + 1)     # quadratic: handled exactly
0.0

julia> distrust(exp(1 ± 0.1))        # gently curved over a narrow spread
0.05029...

julia> distrust(exp(1 ± 1.5))        # wildly non-Gaussian; do not believe it
5.578...
```

This makes `GVar` usable for global optimization, where the mean and radius over a box are used to decide whether the box can contain the optimum. The moments alone will happily miss a narrow, deep minimum — a "slot canyon" — and prune a box that in fact holds the answer. `distrust` gives you a reason to subdivide anyway:

```julia
y = f(lohi(GVar, lo, hi))
if distrust(y) < 1 && mid(y) - 2*rad(y) > best_so_far
    # safe to prune
else
    # bisect regardless of what the parameters say
end
```

`distrust` is a heuristic, not a bound — it cannot make `GVar` rigorous. In particular it cannot see the error caused by treating repeated appearances of a variable as independent, as in `x*x` (write `x^2`). For guaranteed enclosures, use interval arithmetic.

## Relation to other packages

[Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) is a much more fully-developed and featureful package which also handles arithmetic with Gaussian random variables. However, it implements first-order (linear) error propagation, which leads to different rules of arithmetic: compare

```julia
julia> using GaussianRandomVariables

julia> (0 ± 1)^2
1.0 ± 1.4142135623730951 (distrust 0.47)
```

with

```julia
julia> using Measurements

julia> (0 ± 1)^2
0.0 ± 0.0
```

The `GVar` answer is exact: the square of a standard normal is χ²₁, with mean 1 and standard deviation √2. Its distrust is nonzero because χ²₁ is strongly skewed (its skewness is √8), so `mid ∓ rad` is a poor description of where its values actually lie, and the distrust value alerts you to this.

Measurements is recommended for most users, but GaussianRandomVariables can be recommended if second-order accuracy matters in your application.

GaussianRandomVariables is built on top of [ThickNumbers](https://github.com/timholy/ThickNumbers.jl), and the API for working with `GVar`s is described in detail there.
