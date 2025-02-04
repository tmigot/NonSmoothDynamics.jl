struct InterpolationStruct{T1,T2}
  t_vals::Vector{T1}
  x_vals::Matrix{T2}
end

function InterpolationStruct(n::Integer, step_size::T, num_iter::Integer) where {T}
  t_vals = collect(0:step_size:(num_iter*step_size))
  x_vals = Matrix{T}(undef, n, num_iter)
  return InterpolationStruct(t_vals, x_vals)
end

function InterpolationStruct(
  n::Integer,
  t0::AbstractFloat,
  tf::AbstractFloat,
  N::Integer = 100,
)
  step_size = (tf - t0) / (N - 1)
  T = typeof(step_size)
  t_vals = collect(t0:step_size:tf)
  x_vals = Matrix{T}(undef, n, N)
  return InterpolationStruct(t_vals, x_vals)
end
