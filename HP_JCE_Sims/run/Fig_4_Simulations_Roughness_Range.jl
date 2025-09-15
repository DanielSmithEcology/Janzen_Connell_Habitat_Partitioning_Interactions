using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include("../src/HP_JCE_Sims.jl")
using .HP_JCE_Sims
using .HP_JCE_Sims.Simulation
using .HP_JCE_Sims.Parameters
using .HP_JCE_Sims.Environment: match_distribution, generate_habitat
using .HP_JCE_Sims.Simulation: get_neighbors, scaled_JCE_habitat_covariance

using CSV, DataFrames
using Base.Threads: @spawn, @sync
using Random

# === Output directory ===
output_dir = joinpath(@__DIR__, "..", "Outputs", "Fig_4", "simulation_outputs_Fig_4", "Fig_4AB")
isdir(output_dir) || mkpath(output_dir)

# === Parameter values to vary ===
λ_vals = [0.0, 0.1, 0.5, 1.0, 10.0, 100.0]
range_vals = [0.025, 0.05, 0.1, 0.15, 0.2, 0.5]

# === Shared results container ===
results = Vector{NamedTuple{(:species_richness, :Sigma_str, :range_str, :scaled_cov), Tuple{Int, String, String, Float64}}}()
result_lock = Threads.SpinLock()

# === Launch threaded simulations ===
tasks = []

for λ in λ_vals, range_val in range_vals
    push!(tasks, @spawn begin
        config = SimulationConfig(
            grid_size = 175,
            n_species = 500,
            n_generations = 50000,
            sigma = 0.02,
            JC_strength = 0.5,
            dispersal_on = false,
            dispersal_limit = 6,
            mortality_rate = 0.2,
            immigration_rate = 0.3*1e-5,
            dispersal_scale = 5.0,  # required but not used
            Moore_Neigh = 3,
            spatial_autocorr = true,
            roughness_lambda = λ,
            range_param_sim = range_val,
            fecundity_meanlog = 0.0,
            fecundity_sdlog = 0.3,
            temporal_rho = 0.0
        )

        comm, props = run_simulation(config)

        # Run scaled covariance analysis
        optima = collect(range(0, 1, length=config.n_species))
        habitat = generate_habitat(config.grid_size; sill=1.0, range_param=range_val, nugget=λ, seed=config.seed)

        scaled_cov = scaled_JCE_habitat_covariance(
            comm,
            habitat,
            optima,
            config.sigma,
            config.JC_strength,
            config.Moore_Neigh
        )

        # === Format strings ===
        range_str = replace(string(round(range_val, digits=4)), "." => "_")
        λ_str = replace(string(round(λ, digits=4)), "." => "_")

        # === Save outputs ===
        CSV.write(joinpath(output_dir, "comm_range_$(range_str)_lambda_$(λ_str).csv"), DataFrame(comm, :auto))
        CSV.write(joinpath(output_dir, "props_range_$(range_str)_lambda_$(λ_str).csv"), DataFrame(props, :auto))

        # === Store results ===
        richness = length(unique(vec(comm)))
        result = (
            species_richness = richness,
            Sigma_str = λ_str,
            range_str = range_str,
            scaled_cov = scaled_cov
        )

        lock(result_lock) do
            push!(results, result)
        end
    end)
end

# === Wait for all threads to finish ===
@sync begin
    foreach(fetch, tasks)
end

sleep(1)  # Give time for I/O to complete


# === Save final results ===
df = DataFrame(results)
CSV.write(joinpath(output_dir, "Summary_Fig_4AB.csv"), df)



#using Pkg
#Pkg.activate(".")  # activates current directory as project
#cd("C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims")
#include("run/Fig_4_Simulations_Roughness_Range.jl")


