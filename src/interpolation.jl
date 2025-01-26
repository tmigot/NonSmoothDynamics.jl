using Interpolations

# Time vector and discrete solutions
t_vals = collect(0:step_size:(num_iter * step_size))  # Time points
x_vals_t = eachcol(x_vals)  # Convert x_vals (matrix) to a vector of columns

# Create a linear interpolation function for each component of x
interpolations = [LinearInterpolation(t_vals, [x[i] for x in x_vals_t]) for i in 1:length(x0)]

# Interpolated solution
function interpolated_solution(t)
    return [interp(t) for interp in interpolations]
end
