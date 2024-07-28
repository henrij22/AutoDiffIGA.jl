struct BernoulliBeam{𝔅}
    coords::Vector{Float64}
    basis::𝔅

    fext::Function
    EI::Float64

    p::Integer
    dofs::Integer
end

function BernoulliBeam(
        coords::Vector{Float64}, basis::𝔅, f::Function, EI::AbstractFloat) where {𝔅}
    p = basis.degree
    dofs = basis.n

    BernoulliBeam{𝔅}(coords, basis, f, EI, p, dofs)
end

function displacement(beam::BernoulliBeam, d)
    u(ξ) = eval(beam.basis, d, ξ, beam.dofs)
end

function jacobian(beam::BernoulliBeam, d)
    X(ξ) = eval(beam.basis, d, ξ, beam.dofs)
    J(ξ) = ForwardDiff.derivative(X, ξ)
end

function strain(beam::BernoulliBeam, d)
    u = displacement(beam, d)
    J = jacobian(beam, beam.coords)
    Jinv(ξ) = inv(J(ξ))

    u′(ξ) = ForwardDiff.derivative(u, ξ)
    u′′(ξ) = ForwardDiff.derivative(u′, ξ)

    κ(ξ) = u′′(ξ) * Jinv(ξ) * Jinv(ξ)
end

function strainEnergy(beam::BernoulliBeam, d)
    κ = strain(beam, d)
    J = jacobian(beam, beam.coords)
    detJ(ξ) = abs(J(ξ))

    Ψ(ξ) = 0.5 * κ(ξ) * beam.EI * κ(ξ) * detJ(ξ)

    x, w = makeIntegrationRule(numIntegrationPoints(beam.p), unique(beam.basis.knotVec))
    return dot(w, Ψ.(x))
end

function load(beam::BernoulliBeam, d)
    u = displacement(beam, d)

    f(ξ) = u(ξ) * beam.fext(ξ)
end

function loadEnergy(beam::BernoulliBeam, d)
    J = jacobian(beam, d)
    f = load(beam, beam.coords)

    Ψ(ξ) = -f(ξ) * J(ξ)

    x, w = makeIntegrationRule(numIntegrationPoints(beam.p), unique(beam.basis.knotVec))
    return dot(w, Ψ.(x))
end

function computeK(beam::BernoulliBeam)
    x = zeros(beam.dofs)
    ForwardDiff.hessian(d -> strainEnergy(beam, d), x)::Matrix{Float64}
end

function computeF(beam::BernoulliBeam)
    x = zeros(beam.dofs)
    ForwardDiff.gradient(d -> loadEnergy(beam, d), x)::Vector{Float64}
end
