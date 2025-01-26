"""
    project_hyperplane(x::Vector{Float64}, a::Vector{Float64}, b::Float64) -> Vector{Float64}

Projects a vector `x` onto the hyperplane defined by `a^T * x = b`.

# Arguments
- `x::Vector{Float64}`: The vector to project.
- `a::Vector{Float64}`: The normal vector defining the hyperplane.
- `b::Float64`: The offset of the hyperplane.

# Returns
- `Vector{Float64}`: The projected vector.
"""
function project_hyperplane(x::Vector{Float64}, a::Vector{Float64}, b::Float64)
    return x - ((dot(a, x) - b) / dot(a, a)) * a
end

"""
    project_ball(x::Vector{Float64}, c::Vector{Float64}, r::Float64) -> Vector{Float64}

Projects a vector `x` onto the ball centered at `c` with radius `r`.

# Arguments
- `x::Vector{Float64}`: The vector to project.
- `c::Vector{Float64}`: The center of the ball.
- `r::Float64`: The radius of the ball.

# Returns
- `Vector{Float64}`: The projected vector.
"""
function project_ball(x::Vector{Float64}, c::Vector{Float64}, r::Float64)
    dist = norm(x - c)
    return dist <= r ? x : c + r * (x - c) / dist
end

"""
    project_box(x::Vector{Float64}, l::Vector{Float64}, u::Vector{Float64}) -> Vector{Float64}

Projects a vector `x` onto the box defined by lower bounds `l` and upper bounds `u`.

# Arguments
- `x::Vector{Float64}`: The vector to project.
- `l::Vector{Float64}`: The lower bounds of the box.
- `u::Vector{Float64}`: The upper bounds of the box.

# Returns
- `Vector{Float64}`: The projected vector.
"""
function project_box(x::Vector{Float64}, l::Vector{Float64}, u::Vector{Float64})
    return clamp.(x, l, u)
end

"""
    project_simplex(x::Vector{Float64}) -> Vector{Float64}

Projects a vector `x` onto the probability simplex, defined as:
    {y ∈ R^n : y_i ≥ 0, ∑ y_i = 1}

This implementation avoids sorting by iteratively adjusting the elements of `x`
to meet the simplex constraints.

# Arguments
- `x::Vector{Float64}`: The input vector to be projected.

# Returns
- `Vector{Float64}`: The projected vector lying on the simplex.
"""
function project_simplex(x::Vector{Float64})
    n = length(x)
    # Shift x by its average to simplify constraints
    offset = (sum(x) - 1) / n
    x_shifted = x .- offset
    
    # Iteratively adjust negative values to zero while maintaining the constraint
    while any(x_shifted .< 0)
        negative_mask = x_shifted .< 0
        nonnegative_count = sum(.!negative_mask)
        total_negative = sum(x_shifted[negative_mask])
        
        # Distribute the deficit among the nonnegative entries
        offset = total_negative / nonnegative_count
        x_shifted[.!negative_mask] .+= offset
        x_shifted[negative_mask] .= 0
    end
    
    return x_shifted
end

"""
    project_l2_ball(x::Vector{Float64}, r::Float64) -> Vector{Float64}

Projects a vector `x` onto the L2 norm ball of radius `r`.

# Arguments
- `x::Vector{Float64}`: The vector to project.
- `r::Float64`: The radius of the L2 norm ball.

# Returns
- `Vector{Float64}`: The projected vector.
"""
function project_l2_ball(x::Vector{Float64}, r::Float64)
    norm_x = norm(x)
    return norm_x <= r ? x : (r / norm_x) * x
end

"""
    project_positive_orthant(x::Vector{Float64}) -> Vector{Float64}

Projects a vector `x` onto the positive orthant (all elements >= 0).

# Arguments
- `x::Vector{Float64}`: The vector to project.

# Returns
- `Vector{Float64}`: The projected vector.
"""
function project_positive_orthant(x::Vector{Float64})
    return max.(x, 0)
end

"""
    project_halfspace(x::Vector{Float64}, a::Vector{Float64}, b::Float64) -> Vector{Float64}

Projects a vector `x` onto the half-space defined by `a^T * x <= b`.

# Arguments
- `x::Vector{Float64}`: The vector to project.
- `a::Vector{Float64}`: The normal vector defining the half-space.
- `b::Float64`: The offset of the half-space.

# Returns
- `Vector{Float64}`: The projected vector.
"""
function project_halfspace(x::Vector{Float64}, a::Vector{Float64}, b::Float64)
    return dot(a, x) <= b ? x : x - ((dot(a, x) - b) / dot(a, a)) * a
end
