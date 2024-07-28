function map_integration_rules(points, weights, knots::Vector{Float64})
    n = length(points) # Number of integration points
    m = length(knots) - 1 # Number of knot spans
    transformed_points = Vector{Float64}() # To store transformed integration points
    scaled_weights = Vector{Float64}() # To store scaled weights

    for i in 1:m
        a = knots[i]
        b = knots[i + 1]
        mid = (a + b) / 2.0
        half_length = (b - a) / 2.0

        for j in 1:n
            new_point = mid + half_length * points[j]
            new_weight = half_length * weights[j]
            push!(transformed_points, new_point)
            push!(scaled_weights, new_weight)
        end
    end

    return transformed_points, scaled_weights
end

# Define the function to map integration rules to 2D
function map_integration_rules_2d(points, weights, u_knots, v_knots)
    n = length(points) # Number of integration points in 1D
    m_u = length(u_knots) - 1 # Number of knot spans in u direction
    m_v = length(v_knots) - 1 # Number of knot spans in v direction
    transformed_points = Vector{Tuple{Float64, Float64}}() # To store transformed integration points in 2D
    scaled_weights = Vector{Float64}() # To store scaled weights

    for i in 1:m_u
        a_u = u_knots[i]
        b_u = u_knots[i + 1]
        mid_u = (a_u + b_u) / 2.0
        half_length_u = (b_u - a_u) / 2.0

        for k in 1:m_v
            a_v = v_knots[k]
            b_v = v_knots[k + 1]
            mid_v = (a_v + b_v) / 2.0
            half_length_v = (b_v - a_v) / 2.0

            for j in 1:n
                for l in 1:n
                    new_point_u = mid_u + half_length_u * points[j]
                    new_point_v = mid_v + half_length_v * points[l]
                    new_weight = half_length_u * half_length_v * weights[j] * weights[l]
                    push!(transformed_points, (new_point_u, new_point_v))
                    push!(scaled_weights, new_weight)
                end
            end
        end
    end

    return transformed_points, scaled_weights
end

function makeIntegrationRule(k::Integer, uniqueKnots)
    oX, oW = gausslegendre(k)
    map_integration_rules(oX, oW, uniqueKnots)
end

# The knots are assumed to be unique
function makeIntegrationRule(k::Integer, uKnots, vKnots)
    oX, oW = gausslegendre(k)
    map_integration_rules_2d(oX, oW, uKnots, vKnots)
end

function numIntegrationPoints(p)::UInt16
    return ceil((p + 1) / 2)
end
