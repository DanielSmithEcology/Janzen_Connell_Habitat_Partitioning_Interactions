module Simulation



using ..Environment
using ..Competition
using ..SpeciesUtils
using ..Helpers
using ..Parameters: SimulationConfig
import ..Dispersal: dispersal_vector
import ..Competition: recruitment_weights, dispersal_sum_local, dispersal_sum_local_2d, create_2d_kernel, update_dispersal_pressure!, recruitment_weights_tolerance
using Random, Statistics, Distributions, StatsBase, ImageFiltering, NeutralLandscapes



export run_simulation#, scaled_JCE_habitat_covariance, scaled_JCE_habitat_covariance_tolerance, mean_kappa_M


"""
Run the simulation with the given configuration.

Returns:
- Final community matrix (G x G)
- Species proportion matrix (T/n x S) -- sampled every n timesteps of T total timesteps 
"""
# --- in module Simulation ---

function run_simulation(config::SimulationConfig)
    Random.seed!(config.seed)

    G  = config.grid_size
    S  = config.n_species
    T  = config.n_generations
    σ  = config.sigma
    α  = config.JC_strength
    dispersal_on = config.dispersal_on
    r  = config.dispersal_limit
    μ  = config.mortality_rate
    v  = config.immigration_rate
    disp = config.dispersal_scale
    Moore_Nigh = config.Moore_Neigh
    spatial = config.spatial_autocorr
    fecundity_meanlog = config.fecundity_meanlog
    fecundity_sdlog   = config.fecundity_sdlog
    temporal_rho = config.temporal_rho
    λ  = config.roughness_lambda
    range_param_sim = config.range_param_sim

    if spatial
        habitat = Environment.generate_habitat(G; sill=1.0, range_param=range_param_sim, nugget=λ, seed=config.seed)
    else
        habitat = rand(G, G)
    end

    if temporal_rho > 0.0
        size      = (G, G)
        baseline  = rand(G, G)
        habitat_1 = rand(MidpointDisplacement(0.999999999999), size)
        habitat   = match_distribution(habitat_1, baseline)
    end

    # === Species traits ===
    optima         = collect(range(0, 1, length=S))  # unused by tolerance model; harmless
    fecundity_dist = LogNormal(fecundity_meanlog, fecundity_sdlog)

    if config.habitat_model === :tolerance_tradeoff
        fec_draws = rand(fecundity_dist, S)
        fecundity = sort(fec_draws; rev=true)
        ranks = collect(1:S)
    else
        fecundity = rand(fecundity_dist, S)
        ranks = collect(1:S)
    end

    # Dispersal kernel and global fallback (unchanged)
    dispersal_kernel_1d, D_Lim = dispersal_vector(disp, r)
    dispersal_kernel_2d = create_2d_kernel(dispersal_kernel_1d)
    dispersal_pressure  = Array{Float64}(undef, G, G, S)

    # === Initialize state ===
    community = Array{Int}(undef, G, G)
    for i in 1:G, j in 1:G
        hab = habitat[i, j]
        if config.habitat_model === :tolerance_tradeoff
            invS = 1.0 / S
            habitat_effects = Vector{Float64}(undef, S)
            @inbounds for s in 1:S
                habitat_effects[s] = inv(1 + exp(-config.tol_Z * (config.tol_b0 + ranks[s]*invS - hab)))
            end
            ssum = sum(habitat_effects); habitat_effects ./= ssum > 0 ? ssum : 1.0
            community[i, j] = sample(1:S, Weights(habitat_effects))
        else
            habitat_effects = @. exp(-((hab - optima)^2) / (2 * σ^2)) + 1e-3
            habitat_effects ./= sum(habitat_effects)
            community[i, j] = sample(1:S, Weights(habitat_effects))
        end
    end

    # setup for proportions over time (unchanged)
    sample_interval = 25
    n_samples = ceil(Int, T / sample_interval)
    proportions = zeros(n_samples, S)
    sampled_times = Vector{Int}(undef, n_samples)

    # === Run simulation ===
    for t in 1:T
        if temporal_rho > 0.0
            update!(SpatiallyAutocorrelatedUpdater(MidpointDisplacement(.999999999999),0,temporal_rho), habitat)
            habitat = match_distribution(habitat, baseline)
        end

        # Pre-mortality proportions
        p = species_proportions(S, community)

        ### GLOBAL JCE: precompute species-wide global penalty once per gen
        J_global = nothing
        if getfield(config, :JCE_mode) == :global
            M = 2*config.Moore_Neigh + 1
            K = hasproperty(config, :K_total) && config.K_total !== nothing ? config.K_total : (α * M^2)
            J_global = @. exp(-K * p)   # length S
        end

        # Cache pre-mortality community for neighbors / adult counts
        original_community = copy(community)

        if dispersal_on
            update_dispersal_pressure!(dispersal_pressure, community, dispersal_kernel_2d, p, D_Lim)
        end

        # Mortality (unchanged)
        for i in 1:G, j in 1:G
            if rand() < μ
                community[i, j] = 0
            end
        end

        # Recruitment
        for i in 1:G, j in 1:G
            if community[i, j] == 0
                if rand() < v
                    # --- IMMIGRATION ---
                    hab = habitat[i, j]
                    neighbors = get_neighbors(original_community, i, j, Moore_Nigh)

                    if config.habitat_model === :tolerance_tradeoff
                        invS = 1.0 / S
                        habitat_effects = Vector{Float64}(undef, S)
                        @inbounds for s in 1:S
                            habitat_effects[s] = inv(1 + exp(-config.tol_Z * (config.tol_b0 + ranks[s]*invS - hab)))
                        end
                        jc_penalties = if getfield(config, :JCE_mode) == :global
                            J_global
                        else
                            [exp(-α * count(==(s), neighbors)) for s in 1:S]
                        end
                        immigration_weights = fecundity .* habitat_effects .* jc_penalties
                        wsum = sum(immigration_weights); immigration_weights ./= wsum > 0 ? wsum : 1.0
                        community[i, j] = sample(1:S, Weights(immigration_weights))
                        continue
                    else
                        # Gaussian immigration branch
                        habitat_effects = @. exp(-((hab - optima)^2) / (2 * σ^2))
                        jc_penalties = if getfield(config, :JCE_mode) == :global
                            J_global
                        else
                            [exp(-α * count(==(s), neighbors)) for s in 1:S]
                        end
                        immigration_weights = fecundity .* habitat_effects .* jc_penalties
                        wsum = sum(immigration_weights); immigration_weights ./= wsum > 0 ? wsum : 1.0
                        community[i, j] = sample(1:S, Weights(immigration_weights))
                        continue
                    end
                end

                # --- LOCAL RECRUITMENT (lottery) ---
                hab = habitat[i, j]
                neighbors = get_neighbors(original_community, i, j, Moore_Nigh)
                dispersal_weights = dispersal_on ? view(dispersal_pressure, i, j, :)[:] : p

                if config.habitat_model === :tolerance_tradeoff
                    # Under global JCE: remove local JC (α=0) and multiply by J_global
                    if getfield(config, :JCE_mode) == :global
                        weights = Competition.recruitment_weights_tolerance(
                            fecundity, hab, ranks, config.tol_b0, config.tol_Z, neighbors, 0.0;
                            dispersal_weights = dispersal_weights
                        )
                        weights .*= J_global
                    else
                        weights = Competition.recruitment_weights_tolerance(
                            fecundity, hab, ranks, config.tol_b0, config.tol_Z, neighbors, α;
                            dispersal_weights = dispersal_weights
                        )
                    end
                else
                    if getfield(config, :JCE_mode) == :global
                        # Call with α=0 to avoid local JC, then apply global factor
                        weights = Competition.recruitment_weights(
                            fecundity, hab, optima, σ, neighbors, fill(0.0, S);
                            dispersal_weights = dispersal_weights
                        )
                        weights .*= J_global
                    else
                        weights = Competition.recruitment_weights(
                            fecundity, hab, optima, σ, neighbors, fill(α, S);
                            dispersal_weights = dispersal_weights
                        )
                    end
                end

                community[i, j] = sample_species(weights)
            end
        end

        if t % sample_interval == 0
            idx = t ÷ sample_interval
            proportions[idx, :] .= p
            sampled_times[idx] = t
        end
    end

    return community, proportions
end

"""
get_neighbors(community, i, j, r)

Extract the species identities within a square Moore neighborhood of radius `r`
around position (i, j), with wraparound boundaries.
"""
function get_neighbors(community::Matrix{Int}, i::Int, j::Int, r::Int)
    G = size(community, 1)
    neighbors = Int[]
    for dx in -r:r, dy in -r:r
     #   if dx == 0 && dy == 0
     #       continue
     #   end
        ni = mod1(i + dx, G)
        nj = mod1(j + dy, G)
        push!(neighbors, community[ni, nj])
    end
    return neighbors
end

"""
Sample a species from unnormalized recruitment weights.
"""
function sample_species(weights::Vector{Float64})
    probs = weights ./ sum(weights)
    return sample(1:length(probs), Weights(probs))
end

"""
Run multiple simulations with n threads 
"""
function run_multiple_simulations(configs::Vector{SimulationConfig})
    n = length(configs)
    results = Vector{Tuple{Matrix{Int}, Matrix{Float64}}}(undef, n)

    Threads.@threads for i in 1:n
        Random.seed!(configs[i].seed + i)  # unique seed per thread
        results[i] = run_simulation(configs[i])
    end

    return results
end


"""
Compute the mean scaled covariance between JCE strength and habitat suitability
for surviving species in the final community.

Inputs:
- community::Matrix{Int}         : Final G×G community matrix
- habitat::Matrix{Float64}       : G×G matrix of habitat values
- optima::Vector{Float64}        : Species optima (length S)
- σ::Float64                     : Niche width (habitat specialization)
- α::Float64                     : JCE strength
- neigh_radius::Int              : Radius of Moore neighborhood

Returns:
- scaled_covariance::Float64     : Mean covariance (across species) divided by mean spatial fitness
"""
function scaled_JCE_habitat_covariance(
    community::Matrix{Int},
    habitat::Matrix{Float64},
    optima::Vector{Float64},
    σ::Float64,
    α::Float64,
    neigh_radius::Int
)
    G = size(community, 1)
    S = maximum(community)  # Assumes community entries are in 1:S
    surviving_species = unique(community)
    surviving_species = filter(!=(0), surviving_species)  # Remove zeros (empty cells)

    # === A: Compute average spatial fitness ===
    # pho_space_i = ∫ H_i(x) dx  ≈ mean of exp(-((habitat - optima[i])^2) / 2σ^2)
    pho_space_i = [
        mean(@. exp(-((habitat - optima[i])^2) / (2 * σ^2)))
        for i in surviving_species
    ]
    pho_space_avg = mean(pho_space_i)

    # === B: Compute spatial covariance between JCE and habitat effects ===
    covariances = Float64[]
    for s in surviving_species
        jce_effects = Float64[]
        habitat_effects = Float64[]

        for i in 1:G, j in 1:G
            neighbors = get_neighbors(community, i, j, neigh_radius)
            n_neighbors = count(==(s), neighbors)
            jce = 1.0 - exp(-α * n_neighbors)
            hab_effect = exp(-((habitat[i, j] - optima[s])^2) / (2 * σ^2))

            push!(jce_effects, jce)
            push!(habitat_effects, hab_effect)
        end

        cov_val = cov(jce_effects, habitat_effects)
        push!(covariances, cov_val)
    end

    return mean(covariances) / pho_space_avg
end

"""
mean_kappa_M(community, M; threshold=10)

Calculate the average local aggregation index (κ) over all species in the community,
excluding rare species with abundance below `threshold`.

κ_{i,M} = var(n_{i,M}) / mean(n_{i,M}),
where n_{i,M} is the number of individuals of species i in an M×M Moore neighborhood.

Arguments:
- community: G×G matrix of species identities (0 = empty)
- M: neighborhood radius
- threshold: minimum number of individuals to include a species in the calculation

Returns:
- average κ across species with abundance ≥ threshold
"""
function mean_kappa_M(community::Matrix{Int}, M::Int; threshold::Int = 10)
    G = size(community, 1)
    species = unique(vec(community))
    species = filter(x -> x != 0, species)  # remove empty cells (assumed to be 0)

    κ_values = Float64[]

    for sp in species
        # Get all cells where species is present
        presence = community .== sp
        abundance = count(presence)

        if abundance < threshold
            continue
        end

        # For each cell in the grid, count number of sp individuals in M×M neighborhood
        local_counts = Int[]
        for i in 1:G, j in 1:G
            neighbors = get_neighbors(community, i, j, M)
            count_sp = count(==(sp), neighbors)
            push!(local_counts, count_sp)
        end

        μ = mean(local_counts)
        σ² = var(local_counts)
        push!(κ_values, σ² / μ)
    end

    return isempty(κ_values) ? NaN : mean(κ_values)
end



"""
Mean across surviving species of Cov_x(H_k(x)/rho_k, JCE_k(x)) under the
tolerance–fecundity habitat response:

  H_k(x)  = 1 / (1 + exp(-Z * (b0 + r_k/S - x)))
  rho_k   = mean_x H_k(x)
  JCE_k(x)= 1 - exp(-α * n_k(x))

If `ranks` is not provided, we assume species are already ordered by fecundity
as in `:tolerance_tradeoff` and set ranks = 1:S (1 = most fecund).
"""
function scaled_JCE_habitat_covariance_tolerance(
    community::AbstractMatrix{<:Integer},
    habitat::AbstractMatrix{<:Real},
    b0::Float64,
    Z::Float64,
    JC_strength::Float64,
    neigh_radius::Integer;
    ranks::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    jce_mode::Symbol = :pressure,      # :pressure => (1 - exp(-αn)); :survival => exp(-αn)
    rho_eps::Float64 = 1e-12
)::Float64
    G = size(community, 1)
    @assert size(community,2) == G == size(habitat,1) == size(habitat,2)

    surviving_species = filter(!=(0), unique(community))
    isempty(surviving_species) && return NaN

    S = maximum(community)  # assumes species ids are 1..S (0 = empty)
    rks = isnothing(ranks) ? collect(1:S) : ranks
    @assert maximum(surviving_species) <= length(rks) "ranks vector must cover all surviving species"

    covs     = Float64[]
    jce_vec  = Vector{Float64}(undef, G*G)
    Hn_vec   = Vector{Float64}(undef, G*G)

    for s in surviving_species
        rs = rks[s]
        # rho_k = mean_x H_k(x)
        Hbar = mean(@. 1.0 / (1.0 + exp(-Z * (b0 + rs / S - habitat))))
        if Hbar ≤ rho_eps
            continue
        end

        k = 1
        @inbounds for i in 1:G, j in 1:G
            n_s = count(==(s), get_neighbors(community, i, j, neigh_radius))
            jce_vec[k] = (jce_mode === :pressure) ? (1.0 - exp(-JC_strength * n_s)) :
                                                   exp(-JC_strength * n_s)
            x = habitat[i, j]
            H = 1.0 / (1.0 + exp(-Z * (b0 + rs / S - x)))
            Hn_vec[k] = H / Hbar
            k += 1
        end

        push!(covs, cov(Hn_vec, jce_vec))
    end

    return isempty(covs) ? NaN : mean(covs)
end


end # module
