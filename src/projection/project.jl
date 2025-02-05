"""
    project_hyperplane!(sol::Vector{Float64}, x::Vector{Float64}, a::Vector{Float64}, b::Float64) -> Vector{Float64}

Projects a vector `x` onto the hyperplane defined by `a^T * x = b`.

# Arguments
- `sol::Vector{Float64}`: The solution vector to be updated.
- `x::Vector{Float64}`: The vector to project.
- `a::Vector{Float64}`: The normal vector defining the hyperplane.
- `b::Float64`: The offset of the hyperplane.

# Returns
true if the projection is successful.
"""
function project_hyperplane!(sol::S, x::S, a::S, b::T; kwargs...) where {S,T}
  aᵀa = dot(a, a)
  sol .= x - ((dot(a, x) - b) / aᵀa) * a
  return (aᵀa > 0)
end

"""
    project_ball!(sol::Vector{Float64}, x::Vector{Float64}, c::Vector{Float64}, r::Float64) -> Vector{Float64}

Projects a vector `x` onto the ball centered at `c` with radius `r`.

# Arguments
- `sol::Vector{Float64}`: The solution vector to be updated.
- `x::Vector{Float64}`: The vector to project.
- `c::Vector{Float64}`: The center of the ball.
- `r::Float64`: The radius of the ball.

# Returns
true if the projection is successful.
"""
function project_ball!(sol::S, x::S, c::S, r::T; kwargs...) where {S,T}
  dist = norm(x - c)
  sol .= dist <= r ? x : c + r * (x - c) / dist
  return (r >= 0)
end

"""
    project_box!(sol::Vector{Float64}, x::Vector{Float64}, l::Vector{Float64}, u::Vector{Float64}) -> Vector{Float64}

Projects a vector `x` onto the box defined by lower bounds `l` and upper bounds `u`.

# Arguments
- `sol::Vector{Float64}`: The solution vector to be updated.
- `x::Vector{Float64}`: The vector to project.
- `l::Vector{Float64}`: The lower bounds of the box.
- `u::Vector{Float64}`: The upper bounds of the box.

# Returns
true if the projection is successful.
"""
function project_box!(sol::S, x::S, l::S, u::S; kwargs...) where {S}
  sol .= clamp.(x, l, u)
  return all(l .<= u)
end

"""
    project_simplex!(sol::Vector{Float64}, x::Vector{Float64}) -> Vector{Float64}

Projects a vector `x` onto the probability simplex, defined as:
    {y ∈ R^n : y_i ≥ 0, ∑ y_i = 1}

This implementation avoids sorting by iteratively adjusting the elements of `x`
to meet the simplex constraints.

# Arguments
- `sol::Vector{Float64}`: The solution vector to be updated.
- `x::Vector{Float64}`: The input vector to be projected.

# Returns
true if the projection is successful.
"""
function project_simplex!(sol::S, x::S; kwargs...) where {S}
  n = length(x)
  # Shift x by its average to simplify constraints
  offset = (sum(x) - 1) / n
  sol .= x .- offset

  # Iteratively adjust negative values to zero while maintaining the constraint
  while any(sol .< 0)
    negative_mask = sol .< 0
    nonnegative_count = sum(.!negative_mask)
    total_negative = sum(sol[negative_mask])

    # Distribute the deficit among the nonnegative entries
    offset = total_negative / nonnegative_count
    sol[.!negative_mask] .+= offset
    sol[negative_mask] .= 0
  end

  return true
end

"""
    project_l2_ball!(sol::Vector{Float64}, x::Vector{Float64}, r::Float64) -> Vector{Float64}

Projects a vector `x` onto the L2 norm ball of radius `r`.

# Arguments
- `sol::Vector{Float64}`: The solution vector to be updated.
- `x::Vector{Float64}`: The vector to project.
- `r::Float64`: The radius of the L2 norm ball.

# Returns
true if the projection is successful.
"""
function project_l2_ball!(sol::S, x::S, r::T; kwargs...) where {S,T}
  norm_x = norm(x)
  sol .= norm_x <= r ? x : (r / norm_x) * x
  return (r >= 0)
end

"""
    project_positive_orthant!(sol::Vector{Float64}, x::Vector{Float64}) -> Vector{Float64}

Projects a vector `x` onto the positive orthant (all elements >= 0).

# Arguments
- `sol::Vector{Float64}`: The solution vector to be updated.
- `x::Vector{Float64}`: The vector to project.

# Returns
true if the projection is successful.
"""
function project_positive_orthant!(sol::S, x::S; kwargs...) where {S}
  sol .= max.(x, 0)
  return true
end

"""
    project_halfspace!(sol::Vector{Float64}, x::Vector{Float64}, a::Vector{Float64}, b::Float64) -> Vector{Float64}

Projects a vector `x` onto the half-space defined by `a^T * x <= b`.

# Arguments
- `sol::Vector{Float64}`: The solution vector to be updated.
- `x::Vector{Float64}`: The vector to project.
- `a::Vector{Float64}`: The normal vector defining the half-space.
- `b::Float64`: The offset of the half-space.

# Returns
true if the projection is successful.
"""
function project_halfspace!(sol::S, x::S, a::S, b::T; kwargs...) where {S,T}
  aᵀa = dot(a, a)
  sol .= dot(a, x) <= b ? x : x - ((dot(a, x) - b) / aᵀa) * a
  return (aᵀa > 0) | (aᵀa == 0 & all(b .>= 0))
end
