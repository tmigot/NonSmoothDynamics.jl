struct InterpolationStruct{T1,T2}
  t_vals::Vector{T1}
  x_vals::Matrix{T2}
end

function InterpolationStruct(n::Integer, step_size::AbstractFloat, num_iter::Integer)
  t_vals = collect(0:step_size:(num_iter*step_size))
  x_vals = Matrix(undef, n, num_iter)
  return InterpolationStruct(t_vals, x_vals)
end

function InterpolationStruct(n::Integer, t0::T, tf::T, N::Integer = 100) where {T}
  step_size = (tf - t0) / (N - 1)
  t_vals = collect(t0:step_size:tf)
  x_vals = Matrix(undef, n, N)
  return InterpolationStruct(t_vals, x_vals)
end
