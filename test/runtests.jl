module TestMacaulayLab

using Test
using SemialgebraicSets, DynamicPolynomials
using MacaulayLab

function runtests()
    for name in names(@__MODULE__; all = true)
        if !startswith("$(name)", "test_")
            continue
        end
        f = getfield(@__MODULE__, name)
        @testset "$(name)" begin
            f()
        end
    end
end

function test_solve_system()
    # First test from `MacaulayLab/Tests/testsolvesystem.m``
    P = [
        [ 1 2 0 0
         -1 1 1 0
          1 0 0 1],
        [ 2 0 3 0
         -2 1 2 0
         -3 1 1 0],
        [ 1 0 0 3
         -1 1 1 1
         -2 0 0 0],
    ]
    X, output = solvesystem(P)
    @test size(X) == (15, 3)
    @test output["accuracy"] < 1e-10
end

function _in_approx(needle, haystack; kws...)
    return any(x -> isapprox(needle, x; kws...), haystack)
end

function test_solve_system_SemialgebraicSets()
    solver = MacaulayLab.Solver()
    @polyvar x y
    V = @set x^2 + x == 6 && y == x + 1 solver
    sols = collect(V)
    @test length(sols) == 2
    @test _in_approx([2, 3], sols)
    @test _in_approx([-3, -2], sols)
end

end

TestMacaulayLab.runtests()
