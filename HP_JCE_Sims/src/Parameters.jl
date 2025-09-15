module Parameters

export SimulationConfig

using Base: @kwdef
"""
SimulationConfig holds core simulation parameters.

Fields:
- grid_size: Width/height of the square grid (G × G)
- n_species: Number of competing species
- n_generations: Number of simulation steps
- sigma: Niche width for habitat filtering
- JC_strength: Strength of Janzen–Connell competition
- dispersal_on: Toggle for dispersal limitation
- dispersal_limit: Dispersal kernel radius
- mortality_rate: Proportion of adults that die each timestep
- immigration_rate: Probability of regional immigration
- dispersal_scale: Scale parameter for the dispersal kernel
- Moore_Neighborhood: Scale of Janzen-Connell Effects
- spatial_autocorr: Boolean; if false, random habitats; if true, autocorrelated
- kernel_sigma: Standard deviation of Gaussian kernel for smoothing
- roughness_lambda: Blend weight for noise vs. smoothness in habitat
- range_param_sim: Spatial scale of spatial heterogenity
- fecundity_meanlog: mean of lognormal distributed fecundity
- fecundity_sdlog: SD of lognormal distributed fecundity
- temporal_rho: Steps temporal variation; if temporal_rho = 0.0, there's not change over time. 
- seed: Random seed for reproducibility (not to be confused with real seeds) 
- habitat_model: :gaussian (current) | :tolerance_tradeoff
- tol_b0: baseline tolerated proportion (for :tolerance_tradeoff model)
- tol_Z: steepness of the Z-shaped logistic (for :tolerance_tradeoff model
- JCE_mode: determines if JC-effects are local or global (not really JC-effects if global...)
- K_total: can set the K constant for simulations that hold a*M^2 = K; if nothing, computed as JC_strength*(2*Moore_Neigh+1)^2, the same as a* M^2 
- JC_strengths_vec: if not nothing, overrides JC_strength and allows for species-specific JC strengths
"""
# --- in module Parameters ---

Base.@kwdef struct SimulationConfig
    grid_size::Int
    n_species::Int
    n_generations::Int
    sigma::Float64
    JC_strength::Float64
    dispersal_on::Bool
    dispersal_limit::Int
    mortality_rate::Float64
    immigration_rate::Float64
    dispersal_scale::Float64
    Moore_Neigh::Int
    spatial_autocorr::Bool = false
    roughness_lambda::Float64 = 0.1
    range_param_sim::Float64 = 0.15
    fecundity_meanlog::Float64 
    fecundity_sdlog::Float64    
    temporal_rho::Float64 = 0.0  
    seed::Int = 42

    # habitat model selector + params for tolerance model
    habitat_model::Symbol = :gaussian           # :gaussian (current) | :tolerance_tradeoff
    tol_b0::Float64 = 0.05                      # baseline tolerated proportion
    tol_Z::Float64  = 10.0                      # steepness of the Z-shaped logistic

    JCE_mode::Symbol = :local          # :local or :global
    K_total::Union{Nothing,Float64} = nothing  # if nothing, computed as JC_strength*(2*Moore_Neigh+1)^2

    JC_strengths_vec::Union{Nothing,Vector{Float64}} = nothing

end

end
