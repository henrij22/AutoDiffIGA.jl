function uniformRefine(k, p, knots, coordinates)
    for _ in 1:k
        knots, coordinates = refine(knots, coordinates, p, midpoints(unique(knots)))
    end
    knots, coordinates
end

# returns basis, coordinates
function makeGeometry1d(p::Integer, refine::Integer, length::AbstractFloat = 1)
    knots = generateKnotVec(p + 1, p)
    coordinates = [i / p for i in 0:p] .* length

    knots, coordinates = uniformRefine(refine, p, knots, coordinates)
    basis = Bspline(p, knots)
    basis, coordinates
end
