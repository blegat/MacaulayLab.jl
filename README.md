# MacaulayLab

`MacaulayLab.jl` is an interface to the **[MacaulayLab](http://www.macaulaylab.net/)**
solver. It exports the `solvesystem` function that is a thin wrapper on top of the
`solvesystem` MATLAB function and use it to define the `MacaylayLab.Solver` object that
implements the solver-independent `SemialgebraicSets` API for solving algebraic systems.

Consider the system of polynomial equations defined by
`x^2 + x == 6` and `y == x + 1`

To solve it with MacaulayLab directly, do
```julia
julia> using DynamicPolynomial, MacaulayLab

julia> @polyvar x y
(x, y)

julia> P = [x^2 + x - 6, y - x - 1]
2-element Vector{Polynomial{true, Int64}}:
 x² + x - 6
 -x + y - 1

julia> X, output = solvesystem(P);

julia> X
2×2 Matrix{Float64}:
  2.0   3.0
 -3.0  -2.0

julia> output
Dict{String, Any} with 8 entries:
  "residualsbeforeclustering" => [4.66294e-15, 9.76996e-15]
  "accuracybeforeclustering"  => 9.76996e-15
  "time"                      => 1.08337
  "nullities"                 => [0.0, 2.0, 2.0]
  "accuracy"                  => 9.76996e-15
  "residuals"                 => [4.66294e-15, 9.76996e-15]
  "shiftvalues"               => [6.57699, -5.28077]
  "solutions"                 => Any[[6.57699, -5.28077], [2.0, -3.0], [3.0, -2.0]]
```
To solve it with SemialgebraicSets interface, do
```julia
julia> using DynamicPolynomial, MacaulayLab, SemialgebraicSets

julia> solver = MacaulayLab.Solver()
MacaulayLab.Solver()

julia> V = @set x^2 + x == 6 && y == x + 1 solver
Algebraic Set defined by 2 equalities
 x^2 + x - 6.0 = 0
 -x + y - 1.0 = 0

julia> collect(V)
2-element Vector{Vector{Float64}}:
 [2.0000000000000004, 2.999999999999999]
 [-3.0000000000000004, -2.0000000000000004]
```

## Installation

You can install `MacaulayLab.jl` through the Julia package manager:
```julia
] add https://github.com/blegat/MacaulayLab.jl
```
but you first need to make sure that you satisfy the requirements of the
[MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) Julia package and that
the SDPT3 software is installed in your
[MATLAB™](http://www.mathworks.com/products/matlab/) installation.
