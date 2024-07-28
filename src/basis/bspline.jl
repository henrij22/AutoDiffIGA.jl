using NURBS

# Thin wrapper around
struct BSplineBasis
    basis::NURBS.Bspline
    degree::Int64
    knotVec::Vector{Float64}
    n::Int64
end

function BSplineBasis(basis::NURBS.Bspline)
    BSplineBasis(basis, basis.degree, basis.knotVec, numBasisFunctions(basis))
end

function findSpan(b, u, kVec, p)
    if u == kVec[end]
        return b
    end

    # --- Do binary search 
    low = p + 1
    high = b + 1

    mid = Base.midpoint(low, high)

    while u < kVec[mid] || u >= kVec[mid + 1] # is u in interval [ kVec[mid]...kVec[mid+1] ) ?
        if u < kVec[mid]
            high = mid
        else
            low = mid
        end

        mid = Base.midpoint(low, high)
    end

    return mid
end

function basisFun(basis::BSplineBasis, knotSpan, u::T) where {T}
    p = basis.degree
    knotVector = basis.knotVec

    # Allocations
    N = zeros(T, p + 1)
    left = zeros(T, p + 1)
    right = zeros(T, p + 1)

    N[1] = 1.0
    i = knotSpan

    for j in 1:p
        left[j + 1] = u - knotVector[i + 1 - j]
        right[j + 1] = knotVector[i + j] - u
        saved = 0.0

        for r in 1:j
            temp = N[r] / (right[r + 1] + left[j - r + 2])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r + 2] * temp
        end
        N[j + 1] = saved
    end
    return N
end

function eval(basis::BSplineBasis, coefficients, u)
    return eval(basis, coefficients, u, basis.n)
end

function eval(basis::BSplineBasis, coefficients, u, nbasisFun)
    p = basis.degree
    kVec = basis.knotVec
    span = findSpan(nbasisFun, u, kVec, p)

    N = basisFun(basis, span, u)

    result = 0.0

    # span = spans[1]
    for ind in 1:(p + 1)
        result += N[ind] * coefficients[span - p + ind - 1]
    end

    result
end
