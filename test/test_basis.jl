
# Disable Info Logging, Warnings and Errors will still be logged to the console 
disable_logging(Logging.Info)

using NURBS
using LinearAlgebra

function generateBasis(basis, n)
    return [(ξ) -> evalNaive(basis, i, ξ) for i in 1:n]
end
function generateTensorProductBasis(basis, n, m)
    return [((ξ, η) -> evalNaive(basis, j, ξ) * evalNaive(basis, i, η)) for i in 1:n
            for
            j in 1:m]
end

@testset "BSplineBasis" begin
    # test case

    underlyingBasis, _ = makeGeometry1d(2, 0, 1.0)
    basis = BSplineBasis(underlyingBasis)

    coeff = [0.0, 0.3, 0.2]

    N = generateBasis(underlyingBasis, basis.n)

    u(ξ) = dot([N[j](ξ) for j in 1:(basis.n)], coeff)

    testCases = 0:0.033:0.99

    for case in testCases
        @test u(case) ≈ AutoDiffIGA.eval(basis, coeff, case)
    end
end

@testset "BSplineTensorBasis2d" begin
    p = 2
    refineFlag = 1
    basis, _ = makeGeometry2d(p, refineFlag, 1.0, 1.0)

    tensorBasis = BSplineTensorBasis2d(basis, basis)

    dofs = tensorBasis.n
    coeff = collect(range(0, 1, dofs))

    u(ξ) = dot([N[j](ξ...) for j in 1:dofs], coeff)
    N = generateTensorProductBasis(basis, isqrt(dofs), isqrt(dofs))

    testCases = 0:0.033:0.99

    for (i, x) in enumerate(testCases), (j, y) in enumerate(testCases)
        @test u([x, y]) ≈ AutoDiffIGA.eval(tensorBasis, coeff, x, y)
    end
end
