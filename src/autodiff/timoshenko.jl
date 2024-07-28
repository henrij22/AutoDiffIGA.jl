
struct TimoshenkoBeam{𝔅}
    coordsʷ::Vector{Float64}
    basisʷ::𝔅

    coordsᵠ::Vector{Float64}
    basisᵠ::𝔅

    fext::Vector{Function}
    EI::Float64
    GA::Float64

    dofs::Int64
    dofsʷ::Int64
    dofsᵠ::Int64

    indicesʷ::StepRange{Int64, Int64}
    indicesᵠ::StepRange{Int64, Int64}
end

function TimoshenkoBeam(coordsʷ, basisʷ::𝔅, coordsᵠ, basisᵠ::𝔅,
        f::Vector{Function}, EI::Float64, GA::Float64) where {𝔅}
    dofsʷ = basisʷ.n
    dofsᵠ = basisᵠ.n
    dofs = dofsʷ + dofsᵠ

    indicesʷ = 1:dofsʷ
    indicesᵠ = (dofsʷ + 1):dofs

    TimoshenkoBeam{𝔅}(
        coordsʷ, basisʷ, coordsᵠ, basisᵠ, f, EI, GA, dofs, dofsʷ, dofsᵠ, indicesʷ, indicesᵠ)
end

function displacement(beam::TimoshenkoBeam, d)
    u(ξ) = eval(beam.basisʷ, d, ξ, beam.dofsʷ)
end

function rotation(beam::TimoshenkoBeam, d)
    φ(ξ) = eval(beam.basisᵠ, d, ξ, beam.dofsᵠ)
end

function jacobianW(beam::TimoshenkoBeam, d)
    X(ξ) = eval(beam.basisʷ, d, ξ, beam.dofsʷ)
    J(ξ) = ForwardDiff.derivative(X, ξ)
end

function jacobianφ(beam::TimoshenkoBeam, d)
    X(ξ) = eval(beam.basisᵠ, d, ξ, beam.dofsᵠ)
    J(ξ) = ForwardDiff.derivative(X, ξ)
end

function strain(beam::TimoshenkoBeam, d)
    w = displacement(beam, d[beam.indicesʷ])
    φ = rotation(beam, d[beam.indicesᵠ])

    Jʷ = jacobianW(beam, beam.coordsʷ)
    Jinvʷ(ξ) = inv(Jʷ(ξ))

    # This is actually not needed as the geometry description should be the same for w and φ
    # But it has no real performance hit for AutoDiff as this term is constant
    Jᵠ = jacobianφ(beam, beam.coordsᵠ)
    Jinvᵠ(ξ) = inv(Jᵠ(ξ))

    w′(ξ) = ForwardDiff.derivative(w, ξ)
    φ′(ξ) = ForwardDiff.derivative(φ, ξ)

    κ(ξ) = φ′(ξ) * Jinvᵠ(ξ)
    γ(ξ) = w′(ξ) * Jinvʷ(ξ) + φ(ξ)

    # TODO Return as tuple to avoid allocation (?)
    return [κ, γ]
end

function material(beam::TimoshenkoBeam)
    [beam.EI 0
     0 beam.GA]
end

function strainEnergy(beam::TimoshenkoBeam, d)
    ε = strain(beam, d)
    C = material(beam)
    J = jacobianW(beam, beam.coordsʷ)
    detJ(ξ) = abs(J(ξ))

    function Ψ(ξ)
        εₗ = [e(ξ) for e in ε]
        0.5 * εₗ ⋅ (C * εₗ) * detJ(ξ)
    end
    x, w = makeIntegrationRule(
        numIntegrationPoints(beam.basisʷ.degree), unique(beam.basisʷ.knotVec))
    return dot(w, Ψ.(x))
end

function load(beam::TimoshenkoBeam, d)
    w = displacement(beam, d[beam.indicesʷ])
    φ = rotation(beam, d[beam.indicesᵠ])

    vl = [w, φ]

    vlₗ(ξ) = [v(ξ) for v in vl]
    qₗ(ξ) = [v(ξ) for v in beam.fext]

    f(ξ) = vlₗ(ξ) ⋅ qₗ(ξ)
end

function loadEnergy(beam::TimoshenkoBeam, d)
    J = jacobianW(beam, beam.coordsʷ)
    f = load(beam, d)

    Ψ(ξ) = -f(ξ) * J(ξ)

    x, w = makeIntegrationRule(
        numIntegrationPoints(beam.basisʷ.degree), unique(beam.basisʷ.knotVec))
    return dot(w, Ψ.(x))
end

function computeK(beam::TimoshenkoBeam)
    x = zeros(beam.dofs)
    ForwardDiff.hessian(d -> strainEnergy(beam, d), x)::Matrix{Float64}
end

function computeF(beam::TimoshenkoBeam)
    x = zeros(beam.dofs)
    ForwardDiff.gradient(d -> loadEnergy(beam, d), x)
end
