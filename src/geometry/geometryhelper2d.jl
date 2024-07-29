function generateCoordinates(p, length, height)
    p1 = [i / p for i in 0:p] .* height
    p2 = [i / p for i in 0:p] .* length

    coordinates = Vector{Float64}()
    for i in p1, j in p2
        push!(coordinates, j)
        push!(coordinates, i)
    end
    return coordinates
end

function transformToControlPoints(coordinates::AbstractArray, dim)
    controlPoints = Array{SVector{3, Float64}}(undef, dim, dim)

    for 𝓲 in 1:dim, 𝓳 in 1:dim
        k = (𝓲 - 1) * dim + 𝓳
        i = 2 * k - 1
        j = 2 * k

        vec = @SVector [coordinates[i], coordinates[j], 0]
        controlPoints[𝓲, 𝓳] = vec
    end
    controlPoints
end

function transformFromControlPoints(
        controlPoints::AbstractArray{SVector{3, Float64}}, dim::Int)
    num_points = dim * dim
    coordinates = Vector{Float64}(undef, num_points * 2)

    for 𝓲 in 1:dim, 𝓳 in 1:dim
        k = (𝓲 - 1) * dim + 𝓳
        i = 2 * k - 1
        j = 2 * k

        vec = controlPoints[𝓲, 𝓳]
        coordinates[i] = vec[1]
        coordinates[j] = vec[2]
    end

    return coordinates
end

# returns basis, coordinates
function makeGeometry2d(
        p::Integer, refFlag::Integer, length::AbstractFloat = 1, height::AbstractFloat = 1)
    knots = NURBS.generateKnotVec(p + 1, p)
    coordinates = generateCoordinates(p, length, height)

    basis = Bspline(p, knots)
    patch = BsplineSurface(basis, basis, transformToControlPoints(coordinates, p + 1))

    # now refine
    if refFlag == 1
        patch = refine(patch, U = midpoints(unique(knots)), V = midpoints(unique(knots)))
    elseif refFlag == 2
        patch = refFlag(patch, U = midpoints(unique(knots)), V = midpoints(unique(knots)))
        patch = refFlag(patch, U = midpoints(unique(patch.uBasis.knotVec)),
            V = midpoints(unique(patch.vBasis.knotVec)))
    end

    basis = patch.uBasis
    coordinates = transformFromControlPoints(patch.controlPoints, numBasisFunctions(basis))

    basis, coordinates
end
