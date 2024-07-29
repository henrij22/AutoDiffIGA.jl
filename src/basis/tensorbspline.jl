using NURBS

struct BSplineTensorBasis2d
    basisᵘ::BSplineBasis
    basisᵛ::BSplineBasis
    n::Int64
end

function BSplineTensorBasis2d(u_basisᵘ::NURBS.Bspline, u_basisᵛ::NURBS.Bspline)
    basisᵘ = BSplineBasis(u_basisᵘ)
    basisᵛ = BSplineBasis(u_basisᵛ)

    n = basisᵘ.n * basisᵛ.n

    return BSplineTensorBasis2d(basisᵘ, basisᵛ, n)
end

function BSplineTensorBasis2d(basisᵘ::BSplineBasis, basisᵛ::BSplineBasis)
    n = basisᵘ.n * basisᵛ.n

    return BSplineTensorBasis2d(basisᵘ, basisᵛ, n)
end

# function eval(basis::BSplineTensorBasis2d, coefficients, u, v)
#     pᵘ = basis.basisᵘ.degree
#     pᵛ = basis.basisᵛ.degree

#     knotVecᵘ = basis.basisᵘ.knotVec
#     knotVecᵛ = basis.basisᵛ.knotVec

#     uSpan = findSpan(basis.basisᵘ.n, u, knotVecᵘ, pᵘ)
#     Nu = basisFun(basis.basisᵘ, uSpan, u)

#     vSpan = findSpan(basis.basisᵛ.n, v, knotVecᵛ, pᵛ)
#     Nv = basisFun(basis.basisᵛ, vSpan, v)

#     result = 0.0

#     # Coordinates nned to be in a matrix of 
#     coefficients

#     uind = uSpan - pᵘ - 1
#     for i in 1:(pᵛ+1)
#         tmp = 0.0
#         for k in 1:(pᵘ+1)
#             tmp += Nu[k] * coefficients[uind+k]
#         end
#         result += Nv[i] * tmp
#     end
#     result
# end

function eval(basis::BSplineTensorBasis2d, coefficients, u, v)
    pᵘ = basis.basisᵘ.degree
    pᵛ = basis.basisᵛ.degree

    knotVecᵘ = basis.basisᵘ.knotVec
    knotVecᵛ = basis.basisᵛ.knotVec

    uSpan = findSpan(basis.basisᵘ.n, u, knotVecᵘ, pᵘ)
    Nu = basisFun(basis.basisᵘ, uSpan, u)

    vSpan = findSpan(basis.basisᵛ.n, v, knotVecᵛ, pᵛ)
    Nv = basisFun(basis.basisᵛ, vSpan, v)

    result = 0.0

    # Coordinates nned to be in a matrix of 
    coefficients = reshape(coefficients, basis.basisᵘ.n, basis.basisᵛ.n)

    uind = uSpan - pᵘ - 1
    for i in 1:(pᵛ + 1)
        tmp = 0.0
        vind = vSpan - pᵛ + i - 1
        for k in 1:(pᵘ + 1)
            tmp += Nu[k] * coefficients[uind + k, vind]
        end
        result += Nv[i] * tmp
    end
    result
end

# using LinearAlgebra
# p = 2
# knots = generateKnotVec(p + 1, p)
# basis = Bspline(p, knots)

# function generateCoordinates(p, length, height)
#     p1 =  [i / p for i in 0:p] .* height
#     p2 =  [i / p for i in 0:p] .* length

#     coordinates = Vector{Float64}()
#     for i in p1, j in p2
#         push!(coordinates, j)
#         push!(coordinates, i)
#     end
#     return coordinates
# end

# function generateTensorProductBasis(basis, n, m)
#     [((ξ, η) -> evalNaive(basis, j, ξ) * evalNaive(basis, i, η)) for i in 1:n for j in 1:m]
# end

# coordinates = generateCoordinates(p,1, 1)

# n, m = numBasisFunctions(basis), numBasisFunctions(basis)

# dofs_u = n^2
# dofs_v = m^2
# indices_u = 1:2:(dofs_u+dofs_v-1)
# indices_v = 2:2:(dofs_u+dofs_v)

# N = generateTensorProductBasis(basis, n, m)

# jk(ξ) = dot([N[j](ξ...) for j in 1:dofs_u], coordinates[indices_u])

# jk2(ξ) = dot([N[j](ξ...) for j in 1:dofs_u], coordinates[indices_v])

# jk([0.334, 0.4])
# jk2([0.334, 0.4])

# tensorBasis = BSplineTensorBasis2d(basis, basis)

# using StaticArrays

# eval(tensorBasis, coordinates[indices_u], 0.334, 0.4)
# eval(tensorBasis, coordinates[indices_v], 0.334, 0.4)

# mat = zeros(Bool, 11, 11)
# mat2 = zeros(Bool, 11, 11)

# for (i, x) in enumerate(0:0.1:1), (j, y) in enumerate(0:0.1:1)
#     mat[i,j] = eval(tensorBasis, coordinates[indices_u], x, y) ≈ jk([x,y])    
#     mat2[i,j] = eval(tensorBasis, coordinates[indices_v], x, y) ≈ jk2([x,y])    
# end

# @show mat
# @show mat2
