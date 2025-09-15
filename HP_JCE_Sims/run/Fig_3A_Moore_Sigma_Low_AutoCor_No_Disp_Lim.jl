using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include("../src/HP_JCE_Sims.jl")
using .HP_JCE_Sims
using .HP_JCE_Sims.Simulation
using .HP_JCE_Sims.Parameters

using CSV, DataFrames
using Base.Threads: @spawn, nthreads, threadid
using Base.Threads: @threads, @sync
using Base: sleep

# === Output directory ===
output_dir = joinpath(@__DIR__, "..", "Outputs", "Fig_3", "simulation_outputs_Fig_3A")
isdir(output_dir) || mkpath(output_dir)

# === Parameter values to vary ===
Moore_Neigh_vals = [3, 2, 1, 0, "none"]
sigma_vals = [0.02, 0.05, 0.1, 0.2, 9999999]
λ = 10.0  # Fixed roughness value

# === Shared richness results container ===
richness_results = Vector{NamedTuple{(:species_richness, :sigma_str, :lambda_str, :neigh_str), Tuple{Int, String, String, String}}}()
result_lock = Threads.SpinLock()  # or ReentrantLock() would work 


# === Launch threaded simulations ===
tasks = []

for neigh in Moore_Neigh_vals, σ in sigma_vals
    push!(tasks, @spawn begin
        # === Handle JC logic ===
        JanCon_a, neigh_val = neigh == "none" ? (0.0, 0) : (0.5, neigh) # if no JCEs, set a =0.0; otherwise set to 0.5

        config = SimulationConfig(
            grid_size = 175,
            n_species = 500,
            n_generations = 50000, #technically, number of timesteps, not generations per se;  n_generations*mortality_rate ≈ generations
            sigma = σ,
            JC_strength = JanCon_a,
            dispersal_on = false,
            dispersal_limit = 6,
            mortality_rate = 0.2,
            immigration_rate = .3*1e-5,
            dispersal_scale = 5.0,
            Moore_Neigh = neigh_val,
            spatial_autocorr = true,
            roughness_lambda = λ,
            range_param_sim = 0.2,
            fecundity_meanlog = 0.0,
            fecundity_sdlog = 0.3,
            temporal_rho = 0.0
        )

        comm, props = run_simulation(config)

        # === Format strings ===
        neigh_str = string(neigh)
        sigma_str = replace(string(round(σ, digits=4)), "." => "_")
        lambda_str = replace(string(round(λ, digits=4)), "." => "_")

        # === Save main outputs ===
        CSV.write(joinpath(output_dir, "comm_neigh_$(neigh_str)_sigma_$(sigma_str)_lambda_$(lambda_str)_No_Disp_Lim.csv"), DataFrame(comm, :auto))
        CSV.write(joinpath(output_dir, "props_neigh_$(neigh_str)_sigma_$(sigma_str)_lambda_$(lambda_str)_No_Disp_Lim.csv"), DataFrame(props, :auto))

        # === Compute and store richness summary ===
        species_richness = length(unique(vec(comm)))
        result = (
            species_richness = species_richness,
            sigma_str = sigma_str,
            lambda_str = lambda_str,
            neigh_str = neigh_str
        )

        lock(result_lock) do
            push!(richness_results, result)
        end

    end)
end

# === Wait for all threads to finish ===
#foreach(fetch, tasks)
@sync begin
    foreach(fetch, tasks)
end

sleep(1)  # Give time for I/O to complete

# === Save final richness summary ===
df_richness = DataFrame(richness_results)

CSV.write(joinpath(output_dir, "Richness_M_sigmah_No_Disp_Lim_Low_Autocorr.csv"), df_richness)



#using Pkg
#Pkg.activate(".")  # activates current directory as project
#cd("C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims")
#include("run/Fig_3A_Moore_Sigma_Low_AutoCor_No_Disp_Lim.jl")
