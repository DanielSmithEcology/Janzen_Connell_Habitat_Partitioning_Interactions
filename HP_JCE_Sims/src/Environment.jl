module Environment
using ImageFiltering, Random, Statistics, GaussianRandomFields


export habitat_effect, janzen_connell, closest_greater_index, generate_habitat_with_noise, match_distributio, rescale_to_uniform, generate_habitat



"""
    habitat_effect(habitat::Real, sigma::Real, opt::Real)

Compute the Gaussian fitness effect of an individual in a habitat,
given niche width `sigma` and optimum `opt`.
"""
function habitat_effect(habitat::Real, sigma::Real, opt::Real)
    exponent = -((habitat - opt)^2) / (2 * sigma^2)
    return exp(exponent)
end

"""
    janzen_connell(focal::Int, neighborhood::Vector{Int}, a::Vector{Float64})

Compute Janzen–Connell effect: probability of survival is reduced
based on local abundance of conspecifics, scaled by parameter `a[focal]`.
"""
function janzen_connell(focal::Int, neighborhood::Vector{Int}, a::Vector{Float64})
    n_focal = count(x -> x == focal, neighborhood)
    return exp(-a[focal] * n_focal)
end



"""
    generate_habitat_with_noise(G::Int; σ::Float64=1.0, λ::Float64=0.5)

Generate a G×G habitat grid with spatial autocorrelation and adjustable roughness.
- `σ`: standard deviation of Gaussian kernel (controls smoothness).
- `λ`: noise weight (0 = smooth, 1 = rough).

Returns:
- A matrix of Float64 values in [0, 1], matching the distribution of a uniform random field.
"""
function generate_habitat_with_noise(G::Int; σ::Float64=5.0, λ::Float64=0.1)
    base = rand(G, G)
    kernel = Kernel.gaussian((σ, σ))
    smoothed = imfilter(base, kernel, Pad(:circular))
    blended = λ * base .+ (1 - λ) * smoothed

    sorted_blended = sort(vec(blended))
    sorted_base = sort(vec(base))
    flat_blended = reshape(blended, :, 1)
    ranks = [searchsortedfirst(sorted_blended, x) for x in flat_blended]
    matched = reshape(sorted_base[ranks], G, G)

    return matched
end


function match_distribution(source::Matrix{<:Real}, target::Matrix{<:Real})
    flat_source = vec(source)
    flat_target = vec(target)

    sorted_source = sort(flat_source)
    sorted_target = sort(flat_target)

    ranks = sortperm(flat_source)
    matched = similar(flat_source)
    matched[ranks] = sorted_target

    return reshape(matched, size(source))
end


using GaussianRandomFields, Random

"""
    generate_habitat(G; sill=1.0, range=0.2, nugget=0.0, seed=nothing)

Generates a `G x G` habitat matrix as a Gaussian random field with optional nugget effect
and rescales the field to match a uniform distribution.

# Arguments
- `G`: grid size (e.g., 175)
- `sill`: controls the large-scale variance (default 1.0)
- `range`: spatial range parameter `r` 
- `nugget`: nugget effect variance (added as white noise)
- `seed`: random seed (optional, for reproducibility)

# Returns
- `G x G` matrix with values approximately ∼ Uniform(0, 1)
"""
function generate_habitat(G; sill=1.0, range_param=0.2, nugget=0.0, seed=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    # Convert paper's r to Gaussian kernel scale (λ)
    λ = range_param / sqrt(3)

    # 1D coordinates from 0 to 1
    x = range(0, 1, length=G)
    y = range(0, 1, length=G)

    # Covariance function and GRF model
    covfun = CovarianceFunction(2, Gaussian(λ, σ=sqrt(sill)))
    grf = GaussianRandomField(covfun, CirculantEmbedding(), x, y)
    field = sample(grf)

    # Add nugget effect if any
    if nugget > 0.0
        field .+= sqrt(nugget) .* randn(size(field))
    end

    # Rescale to match uniform distribution
    return rescale_to_uniform(field)
end

# Helper function: rescale to uniform(0,1) while preserving rank order
function rescale_to_uniform(field::Matrix{<:Real})
    flat = vec(field)
    sorted = sort(flat)
    ranks = sortperm(flat)
    uniform_target = collect(range(0.0, 1.0, length=length(flat)))
    matched = similar(flat)
    matched[ranks] = uniform_target
    return reshape(matched, size(field))
end



end  # module Environment
