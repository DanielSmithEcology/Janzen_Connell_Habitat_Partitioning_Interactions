using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include("../src/HP_JCE_Sims.jl")
using .HP_JCE_Sims
using .HP_JCE_Sims.Simulation
using .HP_JCE_Sims.Parameters
using .HP_JCE_Sims.Environment: generate_habitat
using .HP_JCE_Sims.Simulation: scaled_JCE_habitat_covariance

using CSV, DataFrames
using Base.Threads: @spawn, @sync
using Random

# ========= Output dir =========
output_dir = joinpath(@__DIR__, "..", "Outputs", "Fig_6", "simulation_outputs_Fig_6", "Fig_6_JCE_HP")
isdir(output_dir) || mkpath(output_dir)

# ========= Experiment design =========
# Hold total JCE-strength K = a*M^2 fixed. You can pass multiple K if you want.
K_vals = [25.0]                     # e.g., try [4.0, 9.0, 16.0] later if you want a robustness sweep
M_vals = [1, 3, 5, 7, 9,11,13,15,17,19,21]           # must be odd (code checks this)

# Optional: habitat structure (set singletons to "hold habitat fixed")
λ_vals     = [0.0]                 # roughness / nugget
range_vals = [0.15]                # spatial range (0-1 scale for your generator)
# If you want a grid, e.g. λ_vals=[0.0,0.5], range_vals=[0.05,0.2], it will do all combos

# ========= Shared containers =========
results = Vector{NamedTuple{(:species_richness, :M, :a, :K, :lambda_str, :range_str, :scaled_cov), 
                            Tuple{Int, Int, Float64, Float64, String, String, Float64}}}()
result_lock = Threads.SpinLock()

# ========= Fixed sim params (tune as you see fit) =========
const GRID          = 175
const N_SPECIES     = 500
const N_GENS        = 50_000
const SIGMA_H       = 0.02
const MORTALITY     = 0.2
const IMM_RATE      = 0.3e-5
const DISP_ON       = false
const DISP_LIM      = 6
const DISP_SCALE    = 5.0          # required but unused if dispersal off
const TEMP_RHO      = 0.0
const SPATIAL_ON    = true

# ========= Threaded run =========
tasks = []

ctr = 0
for K in K_vals, M in M_vals, λ in λ_vals, range_val in range_vals
    if iseven(M) || M < 1
        @warn "Skipping invalid M=$M. It must be odd and ≥1."
        continue
    end
    a = K / (M^2)                # ensure a*M^2 == K
    Moore_Neigh = (M - 1) ÷ 2    # integer division; safe because M is odd
    #ctr += 1

    push!(tasks, @spawn begin
        config = SimulationConfig(
            grid_size          = GRID,
            n_species          = N_SPECIES,
            n_generations      = N_GENS,
            sigma              = SIGMA_H,
            JC_strength        = a,               # <- scaled so that a*M^2 = K
            dispersal_on       = DISP_ON,
            dispersal_limit    = DISP_LIM,
            mortality_rate     = MORTALITY,
            immigration_rate   = IMM_RATE,
            dispersal_scale    = DISP_SCALE,
            Moore_Neigh        = Moore_Neigh,     # <- sets M via 2*Moore_Neigh+1
            spatial_autocorr   = SPATIAL_ON,
            roughness_lambda   = λ,
            range_param_sim    = range_val,
            fecundity_meanlog  = 0.0,
            fecundity_sdlog    = 0.3,
            temporal_rho       = TEMP_RHO
          )

        comm, props = run_simulation(config)

        # Habitat & optima for the scaled covariance analysis
        optima = collect(range(0.0, stop=1.0, length=config.n_species))
        habitat = generate_habitat(config.grid_size; sill=1.0, range_param=range_val, nugget=λ, seed=config.seed)

        
        scaled_cov = scaled_JCE_habitat_covariance(
            comm,
            habitat,
            optima,
            config.sigma,
            config.JC_strength,
            config.Moore_Neigh
        )

        # === Pretty strings ===
        range_str = replace(string(round(range_val, digits=4)), "." => "_")
        λ_str     = replace(string(round(λ,        digits=4)), "." => "_")

        # === Save raw outputs for potential post-hoc analyses ===
        CSV.write(joinpath(output_dir, "comm_K$(K)_M$(M)_range_$(range_str)_lambda_$(λ_str).csv"), DataFrame(comm, :auto))
        CSV.write(joinpath(output_dir, "props_K$(K)_M$(M)_range_$(range_str)_lambda_$(λ_str).csv"), DataFrame(props, :auto))

        # === Store summary ===
        richness = length(unique(vec(comm)))
        result = (
            species_richness = richness,
            M = M,
            a = a,
            K = K,
            lambda_str = λ_str,
            range_str = range_str,
            scaled_cov = scaled_cov
        )

        lock(result_lock) do
            push!(results, result)
        end
    end)
end


@sync begin
    foreach(fetch, tasks)
end

sleep(1)  # Give time for I/O to complete

# ========= Save final summary =========
df = DataFrame(results)
CSV.write(joinpath(output_dir, "Summary_JCE_HP.csv"), df)
println("Saved: ", joinpath(output_dir, "Summary_JCE_HP.csv"))
