module Dispersal

using NumericalIntegration

export dispersal_kernel, dispersal_vector

"""
    dispersal_kernel(x, Az)

Gaussian dispersal kernel with scale parameter `Az`.
Can be broadcast over vectors.
"""
function dispersal_kernel(x::Real, Az::Real)
    return (2 * π * x / (π * Az^2)) * exp(-x^2 / Az^2)
end

"""
    dispersal_vector(Az, maxlen)

Compute a binned dispersal vector from a Gaussian kernel with scale `Az`,
extending up to `maxlen` distance shells.

Returns:
- `disp_vec`: Vector of dispersal probabilities at each distance shell.
- `total_prob`: Total cumulative probability captured.
"""
function dispersal_vector(Az::Float64, maxlen::Int)
    dx = 1e-4
    disp_vec = Float64[]

    # First bin: from 0 to 0.5
    x0 = 0:dx:0.5
    y0 = dispersal_kernel.(x0, Az)
    bin_prob = integrate(x0, y0)
    push!(disp_vec, bin_prob)
    total_prob = bin_prob

    d = 1
    while d < maxlen && total_prob < 0.99
        x = (d - 0.5):dx:(d + 0.5)
        y = dispersal_kernel.(x, Az)
        bin_prob = (1 / (8 * d)) * integrate(x, y)

        push!(disp_vec, bin_prob)
        total_prob += bin_prob*(d*8)
        d += 1
    end

    return disp_vec, total_prob
end

"""
    create_2d_kernel(kernel::Vector{Float64}) -> Matrix{Float64}

Expands a 1D Chebyshev-based dispersal kernel into a full 2D kernel matrix.
Assumes dispersal depends on Chebyshev distance.
"""
function create_2d_kernel(kernel::Vector{Float64})
    r = length(kernel) - 1
    K = zeros(Float64, 2r+1, 2r+1)
    for dx in -r:r, dy in -r:r
        d = max(abs(dx), abs(dy))
        K[dx + r + 1, dy + r + 1] = kernel[d + 1]
    end
    return K
end





end  # module Dispersal
