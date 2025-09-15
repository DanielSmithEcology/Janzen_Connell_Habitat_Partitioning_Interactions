using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include("../src/HP_JCE_Sims.jl")
using .HP_JCE_Sims
using .HP_JCE_Sims.Simulation
using .HP_JCE_Sims.Parameters

using CSV, DataFrames
using Base.Threads: @spawn, nthreads, threadid, @threads, @sync
using Base: sleep

# === Output directory ===
# NOTE: Using the exact folder name you gave ("Llim" with two Ls).
output_dir = joinpath(@__DIR__, "..", "Outputs", "Appendix_SST", "High_Auto_No_Disp_Lim")
isdir(output_dir) || mkpath(output_dir)

# === Parameter values to vary ===
Moore_Neigh_vals = [3, 2, 1, 0, "none"]      # as before (M=0 or "none" disables JCE)
Z_vals           = [0.0, 15.0, 25.0, 50.0, 150.0]
λ = 0.0  # Fixed roughness value

# === Shared richness results container ===
richness_results = Vector{NamedTuple{(:species_richness, :Z_str, :lambda_str, :neigh_str), Tuple{Int, String, String, String}}}()
result_lock = Threads.SpinLock()

# === Launch threaded simulations ===
tasks = []

for neigh in Moore_Neigh_vals, Z in Z_vals
    push!(tasks, @spawn begin
        # === Handle JC logic ===
        JanCon_a, neigh_val = neigh == "none" ? (0.0, 0) : (0.5, neigh)

        config = SimulationConfig(
            grid_size = 175,
            n_species = 500,
            n_generations = 25000,     # generations ≈ n_generations * mortality_rate
            sigma = 0.05,               # unused under :tolerance_tradeoff; kept for compatibility
            JC_strength = JanCon_a,
            dispersal_on = false,
            dispersal_limit = 6,
            mortality_rate = 0.2,
            immigration_rate = 0.3e-5,
            dispersal_scale = 5.0,
            Moore_Neigh = neigh_val,
            spatial_autocorr = true,
            roughness_lambda = λ,
            range_param_sim = 0.5,
            fecundity_meanlog = 0.0,
            fecundity_sdlog = 2.0,
            temporal_rho = 0.0,
            # === NEW: toggle stress-tolerance model ===
            habitat_model = :tolerance_tradeoff,
            tol_b0 = 0.05,
            tol_Z  = Z
        )

        comm, props = run_simulation(config)

        # === Format strings ===
        neigh_str  = string(neigh)
        Z_str      = replace(string(round(Z, digits=4)), "." => "_")
        lambda_str = replace(string(round(λ, digits=4)), "." => "_")

        # === Save main outputs ===
        CSV.write(joinpath(output_dir, "comm_neigh_$(neigh_str)_Z_$(Z_str)_lambda_$(lambda_str)_No_Disp_Lim.csv"),
                  DataFrame(comm, :auto))
        CSV.write(joinpath(output_dir, "props_neigh_$(neigh_str)_Z_$(Z_str)_lambda_$(lambda_str)_No_Disp_Lim.csv"),
                  DataFrame(props, :auto))

        # === Compute and store richness summary ===
        species_richness = length(unique(vec(comm)))  # (kept consistent with your current approach)
        result = (
            species_richness = species_richness,
            Z_str = Z_str,
            lambda_str = lambda_str,
            neigh_str = neigh_str
        )

        lock(result_lock) do
            push!(richness_results, result)
        end
    end)
end

# === Wait for all threads to finish ===
@sync foreach(fetch, tasks)
sleep(1)  # allow I/O to flush

# === Save final richness summary ===
df_richness = DataFrame(richness_results)
CSV.write(joinpath(output_dir, "Richness_M_Z_No_Disp_Lim_High_Autocorr.csv"), df_richness)


#using Pkg
#Pkg.activate(".")  # activates current directory as project
#cd("C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims")
#include("run/Appendix_SST_High_AutoCor_No_Disp_Lim.jl")
