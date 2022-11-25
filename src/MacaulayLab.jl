module MacaulayLab

import MATLAB
import MultivariatePolynomials
const MP = MultivariatePolynomials
import SemialgebraicSets

export solvesystem

function solvesystem(P::Vector{Matrix{Float64}})
    return MATLAB.mxcall(:solvesystem, 2, P)
end

solvesystem(P::Vector{<:Matrix{<:Integer}}) = solvesystem(float.(P))

function _mat(poly::MP.AbstractPolynomialLike)
    P = Matrix{Float64}(undef, MP.nterms(poly), 1 + MP.nvariables(poly))
    for (i, t) in enumerate(MP.terms(poly))
        P[i, 1] = MP.coefficient(t)
        P[i, 2:end] = MP.exponents(t)
    end
    return P
end

function solvesystem(P::Vector{<:MP.AbstractPolynomialLike})
    # hack to make sure they all have all variables:
    P = P .* one(prod(MP.variables(P)))
    return solvesystem(_mat.(P))
end

struct Solver <: SemialgebraicSets.AbstractAlgebraicSolver end

function SemialgebraicSets.solvealgebraicequations(
    V::SemialgebraicSets.AbstractAlgebraicSet,
    ::Solver,
)
    p = SemialgebraicSets.equalities(V)
    X, _ = solvesystem(p)
    return [X[i, :] for i in axes(X, 1)]
end

end # module MacaulayLab
