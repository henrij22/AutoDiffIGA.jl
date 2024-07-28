struct LinearElastic{𝔅}
    coords::Vector{Float64}
    basis::𝔅

    fext::Vector{Function}
    E::Float64
    ν::Float64

    dofs::Int64
    dofsᵘ::Int64
    dofsᵛ::Int64

    indicesᵘ::StepRange{Int64, Int64}
    indicesᵛ::StepRange{Int64, Int64}
end

function LinearElastic(
        coords, basis::𝔅, f::Vector{Function}, E::Float64, ν::Float64) where {𝔅}
    dofs_u = basis.n
    dofs_v = basis.n

    dofs = dofs_u + dofs_v
    indices_u = 1:2:(dofs_u + dofs_v - 1)
    indices_v = 2:2:(dofs_u + dofs_v)

    LinearElastic{𝔅}(coords, basis, f, E, ν, dofs, dofs_u, dofs_v, indices_u, indices_v)
end

function displacement(plane::LinearElastic, d::AbstractArray)
    u(ξ) = eval(plane.basis, d, ξ...)
end

function geometry(plane::LinearElastic, d::AbstractArray)
    x(ξ) = eval(plane.basis, d, ξ...)
end

function jacobian(plane::LinearElastic)
    geo1 = geometry(plane, plane.coords[plane.indicesᵘ])
    geo2 = geometry(plane, plane.coords[plane.indicesᵛ])

    X(ξ) = [geo1(ξ), geo2(ξ)]
    J(ξ) = ForwardDiff.jacobian(X, ξ)
end

function strain(plane::LinearElastic, d::AbstractArray)
    u = displacement(plane, d[plane.indicesᵘ])
    v = displacement(plane, d[plane.indicesᵛ])

    J = jacobian(plane)
    Jinv(ξ) = inv(J(ξ)')

    gradientᵤ(ξ) = ForwardDiff.gradient(u, ξ)
    gradientᵥ(ξ) = ForwardDiff.gradient(v, ξ)

    # Isoparametric concept (geometry transformations)
    # Notation: derivative of u in direction of ξ -> uξ
    uξ(ξ) = gradientᵤ(ξ)[1]
    uη(ξ) = gradientᵤ(ξ)[2]

    uˣ(ξ) = uξ(ξ) * Jinv(ξ)[1, 1] + uη(ξ) * Jinv(ξ)[1, 2]
    uʸ(ξ) = uξ(ξ) * Jinv(ξ)[2, 1] + uη(ξ) * Jinv(ξ)[2, 2]

    vξ(ξ) = gradientᵥ(ξ)[1]
    vη(ξ) = gradientᵥ(ξ)[2]

    vˣ(ξ) = vξ(ξ) * Jinv(ξ)[1, 1] + vη(ξ) * Jinv(ξ)[1, 2]
    vʸ(ξ) = vξ(ξ) * Jinv(ξ)[2, 1] + vη(ξ) * Jinv(ξ)[2, 2]

    εˣ(ξ) = uˣ(ξ)
    εʸ(ξ) = vʸ(ξ)
    εˣʸ(ξ) = uʸ(ξ) + vˣ(ξ)
    return [εˣ, εʸ, εˣʸ]
end

function material(plane::LinearElastic)
    plane.E / (1 - plane.ν^2) * [1 plane.ν 0
                                 plane.ν 1 0
                                 0 0 (1 - plane.ν)/2]
end

function strainEnergy(plane::LinearElastic, d::AbstractArray)
    ε = strain(plane, d)
    J = jacobian(plane)
    C = material(plane)

    # Ψ = ½ ε ⋅ C ⋅ ε
    function Ψ(ξ)
        εₗ = [e([ξ[1], ξ[2]]) for e in ε]
        0.5 * εₗ ⋅ (C * εₗ) * det(J([ξ[1], ξ[2]]))
    end

    x, w = makeIntegrationRule(
        numIntegrationPoints(2p), unique(plane.basis.basisᵘ.knotVec),
        unique(plane.basis.basisᵛ.knotVec))
    return dot(w, Ψ.(x))
end

function volumeLoad(plane::LinearElastic, d::AbstractArray)
    u(ξ) = displacement(plane, d[plane.indicesᵘ])(ξ)
    v(ξ) = displacement(plane, d[plane.indicesᵛ])(ξ)

    vl = [u, v]

    vlₗ(ξ) = [v(ξ) for v in vl]
    qₗ(ξ) = [v(ξ) for v in q]

    vll(ξ) = vlₗ(ξ) ⋅ qₗ(ξ)
end

function volumeLoadEnergy(plane::LinearElastic, d::AbstractArray)
    J = jacobian(plane)
    vl = volumeLoad(plane, d)

    # Ψ = N ⋅ f
    Ψ(ξ) = -vl([ξ[1], ξ[2]]) * det(J([ξ[1], ξ[2]]))

    x, w = makeIntegrationRule(
        numIntegrationPoints(2p), unique(plane.basis.basisᵘ.knotVec),
        unique(plane.basis.basisᵛ.knotVec))
    return dot(w, Ψ.(x))
end

function computeK(plane::LinearElastic)
    x = zeros(plane.dofs)
    ForwardDiff.hessian(d -> strainEnergy(plane, d), x)::Matrix{Float64}
end

function computeF(plane::LinearElastic)
    x = zeros(plane.dofs)
    ForwardDiff.gradient(d -> volumeLoadEnergy(plane, d), x)::Vector{Float64}
end
