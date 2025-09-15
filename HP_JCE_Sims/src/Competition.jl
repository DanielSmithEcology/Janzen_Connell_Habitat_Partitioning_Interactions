module Competition

using ..Helpers
using Statistics, ImageFiltering 

export recruitment_weights, dispersal_sum_local, create_2d_kernel, dispersal_sum_local_2d, update_dispersal_pressure!, recruitment_weights_tolerance


"""
Compute recruitment weights based on species-specific traits, habitat matching, Janzen-Connell effects,
and optionally dispersal weights.

Arguments:
- fecundity: Vector of fecundity values per species
- habitat_value: Local habitat value
- optima: Vector of habitat optima per species
- sigma: Niche width (Gaussian tolerance)
- neighbors: List of species in the local neighborhood
- JC_strengths: Janzen–Connell effect strengths per species
- dispersal_weights: Optional vector of dispersal weights per species

Returns:
- Vector of recruitment weights per species
"""
function recruitment_weights(
    fecundity::Vector{Float64},
    habitat_value::Float64,
    optima::Vector{Float64},
    sigma::Float64,
    neighbors::Vector{Int},
    JC_strengths::Vector{Float64};
    dispersal_weights::Vector{Float64}=ones(length(fecundity))
)
    S = length(fecundity)
    weights = zeros(Float64, S)
    for s in 1:S
        habitat_effect = exp(-((habitat_value - optima[s])^2) / (2 * sigma^2))
        jc_penalty = exp(-JC_strengths[s] * count(==(s), neighbors))
        weights[s] = fecundity[s] * habitat_effect * jc_penalty * dispersal_weights[s]
    end
    return weights
end



"""
Compute recruitment weights under a tolerance–fecundity trade-off.

H_i(x) = 1 / (1 + exp(-Z * (b0 + (rank_i / S) - x)))

Args:
- fecundity::Vector{Float64}      : per-species fecundity (rank 1 is highest)
- habitat_value::Float64          : x in [0,1]
- ranks::Vector{Int}              : length S, rank 1..S (1 = most fecund)
- b0::Float64                     : baseline tolerated proportion
- Z::Float64                      : steepness of the logistic
- neighbors::Vector{Int}          : species in local neighborhood
- JC_strengths                    : Float64 OR Vector{Float64}
- dispersal_weights::Vector{Float64} (kwarg): per-species weights

Returns: Vector{Float64} of recruitment weights (length S)
"""
function recruitment_weights_tolerance(
    fecundity::Vector{Float64},
    habitat_value::Float64,
    ranks::Vector{Int},
    b0::Float64,
    Z::Float64,
    neighbors::Vector{Int},
    JC_strengths;
    dispersal_weights::Vector{Float64}=ones(length(fecundity))
)
    S = length(fecundity)
    @assert length(ranks) == S "ranks must have length S"
    @assert length(dispersal_weights) == S "dispersal_weights must have length S"

    invS = 1.0 / S
    weights = zeros(Float64, S)

    @inbounds for s in 1:S
        αs = (JC_strengths isa AbstractVector) ? JC_strengths[s] : JC_strengths
        # Z-shaped tolerance curve (higher rank tolerates higher x)
        H  = inv(1 + exp(-Z * (b0 + ranks[s]*invS - habitat_value)))
        J  = exp(-αs * count(==(s), neighbors))
        weights[s] = fecundity[s] * H * J * dispersal_weights[s]
    end
    return weights
end



"""
Compute local dispersal contributions + global fallback for each species.

Arguments:
- s: species index
- i, j: focal cell coordinates
- community: G×G matrix of species
- kernel: dispersal kernel (flattened)
- p: global species frequencies
- D_Lim: dispersal limitation strength (0 to 1)

Returns:
- Total weight for species s at (i, j) due to local dispersal + global fallback
"""
function dispersal_sum_local(
    s::Int,
    i::Int,
    j::Int,
    community::Matrix{Int},
    kernel::Vector{Float64},  # dispersal_kernel[d+1]
    p::Vector{Float64},
    D_Lim::Float64
)
    G = size(community, 1)
    len = length(kernel)
    half = len-1  # max Chebyshev distance covered by kernel
    total = 0.0

    for dx in -half:half, dy in -half:half
        ni = mod1(i + dx, G)
        nj = mod1(j + dy, G)

        if community[ni, nj] == s
            d = max(abs(dx), abs(dy))  # Chebyshev distance
            if d + 1 <= length(kernel)
                total += kernel[d + 1]
            end
        end
    end

    return total + (1 - D_Lim) * p[s]  # global fallback
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



function dispersal_sum_local_2d(
    s::Int,
    i::Int,
    j::Int,
    community::Matrix{Int},
    kernel_2d::Matrix{Float64},
    p::Vector{Float64},
    D_Lim::Float64
)
    G = size(community, 1)
    r = (size(kernel_2d, 1) - 1) ÷ 2
    total = 0.0

    for dx in -r:r, dy in -r:r
        ni = mod1(i + dx, G)
        nj = mod1(j + dy, G)

        if community[ni, nj] == s
            weight = kernel_2d[dx + r + 1, dy + r + 1]
            total += weight
        end
    end

    return total + (1 - D_Lim) * p[s]
end


"""
    update_dispersal_pressure!(dispersal_pressure, community, dispersal_kernel_2d, p, D_Lim)

Updates the dispersal pressure at each patch for each species based on local filtering via 2D convolution.
Includes a global fallback proportional to (1 - D_Lim) * p[s].

Arguments:
- `dispersal_pressure::Array{Float64, 3}`: (G×G×S) tensor to be updated in-place.
- `community::Matrix{Int}`: current community matrix of species identities.
- `dispersal_kernel_2d::Matrix{Float64}`: 2D convolution kernel for local dispersal.
- `p::Vector{Float64}`: current species proportions.
- `D_Lim::Float64`: weight of local dispersal vs. global fallback.
"""
function update_dispersal_pressure!(
    dispersal_pressure::Array{Float64, 3},
    community::Matrix{Int},
    dispersal_kernel_2d::Matrix{Float64},
    p::Vector{Float64},
    D_Lim::Float64
)
    G = size(community, 1)
    S = length(p)

    # Identify species still present
    present_species = unique(vec(community))

    for s in present_species
        presence_map = map(x -> x == s ? 1.0 : 0.0, community)
        local_dispersal = imfilter(presence_map, dispersal_kernel_2d, Pad(:circular))
        dispersal_pressure[:, :, s] .= D_Lim * local_dispersal .+ (1 - D_Lim) * p[s]
    end

    # For extinct species: just assign global fallback
    extinct_species = setdiff(1:S, present_species)
    for s in extinct_species
        dispersal_pressure[:, :, s] .= (1 - D_Lim) * p[s]
    end
end




end


