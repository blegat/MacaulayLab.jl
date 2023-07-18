module MacaulayLab

import MATLAB
import MultivariatePolynomials as MP
import SemialgebraicSets as SS

export solvesystem

Base.@kwdef struct Options
    algorithm::String = "null" # selects the algorithm "column", "null" or "cheb"
    blocked::Bool = true # selects the blocked version of the algorithm
    recursive::String = "sparse" # selects an approach to build the null space ("full", "iterative", "recursive", "sparse").
    tol::Float64 = 1e-10 # tolerance in the rank checks
    # TODO shiftpoly # [double(m,n+1)] main shift polynomial (default = [randn(n,1) eye(n)]).
    verbose::Bool = false # displays additional information (default = false).
    clustering::Bool = true # choice to apply clustering to gather the multiplicities
    clustertol::Float64 = 1e-3 # tolerance in the clustering step
    posdim::Bool = false # flag to solve problems with a positive-dimensional solution set at infinity
end

function solvesystem(P::Vector{Matrix{Float64}}, maxdegree::Float64, options::Options = Options())
    params = Dict(f => getfield(options, f) for f in fieldnames(typeof(options)))
    return MATLAB.mxcall(:solvesystem, 2, P, maxdegree, params)
end

function solvesystem(P::Vector{Matrix{Float64}}, maxdegree::Integer, args...)
    return solvesystem(P, Float64(maxdegree), args...)
end

solvesystem(P::Vector{<:Matrix{<:Integer}}, args...) = solvesystem(float.(P), args...)

_vec(t::Tuple) = Float64[t...]
_vec(v::AbstractVector) = v

function _mat(poly::MP.AbstractPolynomialLike)
    P = Matrix{Float64}(undef, MP.nterms(poly), 1 + MP.nvariables(poly))
    for (i, t) in enumerate(MP.terms(poly))
        P[i, 1] = MP.coefficient(t)
        P[i, 2:end] = _vec(MP.exponents(t))
    end
    return P
end

function _mat(P::Vector{<:MP.AbstractPolynomialLike})
    # hack to make sure they all have all variables:
    P = P .* one(prod(MP.variables(P)))
    return _mat.(P)
end

function solvesystem(P::Vector{<:MP.AbstractPolynomialLike}, maxdegree::Real, args...)
    return solvesystem(_mat(P), maxdegree, args...)
end

function solvesystem(P::Vector{<:MP.AbstractPolynomialLike}, ::Nothing, args...)
    return solvesystem(P, args...)
end

function solvesystem(P::Vector{<:MP.AbstractPolynomialLike}, args...)
    return solvesystem(P, sum(MP.maxdegree, P) - length(P) + 2, args...)
end

struct Solver <: SS.AbstractAlgebraicSolver
    maxdegree
    options::Options
end
Solver(maxdegree=nothing) = Solver(maxdegree, Options())

SS.default_gröbner_basis_algorithm(::Any, ::Solver) = SS.NoAlgorithm()

SS.promote_for(::Type{T}, ::Type{Solver}) where {T} = float(T)

function _real(x::Complex{T}) where T
    if abs(imag(x)) > Base.rtoldefault(T)
        @warn("Ignoring imaginary part of `$x`")
    end
    return real(x)
end

_real(x::Real) = x

function _real(x::Vector)
    return _real.(x)
end

function SS.solve(V::SS.AbstractAlgebraicSet, s::Solver)
    X, _ = solvesystem(SS.equalities(V), s.maxdegree, s.options)
    return [_real(X[i, :]) for i in axes(X, 1)]
end

end # module MacaulayLab
