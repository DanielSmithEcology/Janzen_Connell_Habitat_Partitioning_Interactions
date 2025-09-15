using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include("../src/HP_JCE_Sims.jl")
using .HP_JCE_Sims
using .HP_JCE_Sims.Simulation
using .HP_JCE_Sims.Parameters

using CSV, DataFrames
using Base.Threads: @spawn, @sync
using Base: sleep
using Random

# ============== Output dir ==============
output_dir = joinpath(@__DIR__, "..", "Outputs", "Fig_SpeciesPoolSize", "simulation_outputs_Fig_Nspecies")
isdir(output_dir) || mkpath(output_dir)

# ============== Sweep over n_species ==============
#n_species_vals = [15,31, 62, 125, 250, 500, 1000, 2000]

n_species_vals = [10, 65, 125, 250, 500, 750, 1000]


# ============== Build the exact parameter cases you requested ==============
# 1) σ = 9999999 with a sweep of Moore_Neigh in [3,2,1,0,"none"] (λ defaults to reference λ = 0.0)
cases = NamedTuple[]

# 1) σ = 9_999_999 (global) with Moore_Neigh = 3 (λ = 0.0)
push!(cases, (label = "No_HP_sigma_M3", sigma = 9_999_999.0, neigh = 3, lambda = 0.0))

# 2) σ = 0.02 and Moore_Neigh = "none" (λ = 0.0 by reference)
push!(cases, (label = "sigma_0p02_none", sigma = 0.02, neigh = "none", lambda = 0.0))

# 3) σ = 0.02, Moore_Neigh = 3, with λ = 0.0, 0.5, 10.0 (three separate cases)
for lam in (0.0, 0.1, 10.0)
    push!(cases, (label = "sigma_0p02_M3_lam_$(replace(string(lam), "." => "p"))",
                  sigma = 0.02, neigh = 3, lambda = lam))
end

# ============== Helpers for strings ==============
to_str(x) = replace(string(x), "." => "_")
sigma_tag(σ) = (σ ≥ 1e6 ? "No_HP" : to_str(round(σ, digits=6)))
lambda_tag(λ) = to_str(round(λ, digits=6))
neigh_tag(n)  = string(n)

# ============== Shared defaults (from your reference code) ==============
const GRID_SIZE        = 175
const N_GENERATIONS    = 50_000
const DISPERSAL_ON     = false
const DISPERSAL_LIMIT  = 6
const MORTALITY_RATE   = 0.2
const IMMIGRATION_RATE = 0.3e-5 #0.3e-5
const DISPERSAL_SCALE  = 5.0
const SPATIAL_AUTOCORR = true
const RANGE_PARAM_SIM  = 0.5
const FEC_MEANLOG      = 0.0
const FEC_SDLOG        = 0.3
const TEMPORAL_RHO     = 0.0
const JC_DEFAULT       = 0.5   # when JC is "on" (Moore_Neigh is Int)

# ============== Results container ==============
const RowT = NamedTuple{
    (:species_richness, :sigma_str, :lambda_str, :neigh_str, :nsp_str, :case_label),
    Tuple{Int, String, String, String, String, String}
}
richness_results = Vector{RowT}()
result_lock = Threads.SpinLock()

# ============== Parallel runs (single block) ==============
tasks = []

@sync begin
    for nsp in n_species_vals, cs in cases
        push!(tasks, @spawn begin
            # Seed per task for reproducibility
            Random.seed!(hash((nsp, cs.label, cs.sigma, cs.neigh, cs.lambda)))

            # JC logic
            neigh_val = cs.neigh == "none" ? 0 : Int(cs.neigh)
            jc_strength = (cs.neigh == "none") ? 0.0 : JC_DEFAULT

            config = SimulationConfig(
                grid_size = GRID_SIZE,
                n_species = nsp,
                n_generations = N_GENERATIONS,   # ≈ (gens) since gens * mortality_rate ≈ effective gens
                sigma = cs.sigma,
                JC_strength = jc_strength,
                dispersal_on = DISPERSAL_ON,
                dispersal_limit = DISPERSAL_LIMIT,
                mortality_rate = MORTALITY_RATE,
                immigration_rate = IMMIGRATION_RATE,
                dispersal_scale = DISPERSAL_SCALE,
                Moore_Neigh = neigh_val,
                spatial_autocorr = SPATIAL_AUTOCORR,
                roughness_lambda = cs.lambda,
                range_param_sim = RANGE_PARAM_SIM,
                fecundity_meanlog = FEC_MEANLOG,
                fecundity_sdlog = FEC_SDLOG,
                temporal_rho = TEMPORAL_RHO
            )

            comm, props = run_simulation(config)

            # Filenames
            neigh_str   = neigh_tag(cs.neigh)
            sigma_str   = sigma_tag(cs.sigma)
            lambda_str  = lambda_tag(cs.lambda)
            nsp_str     = string(nsp)
            stub = "case_$(cs.label)_nsp_$(nsp_str)_neigh_$(neigh_str)_sigma_$(sigma_str)_lambda_$(lambda_str)"

            # Save outputs
            CSV.write(joinpath(output_dir, "comm_$(stub)_No_Disp_Lim.csv"),  DataFrame(comm, :auto))
            CSV.write(joinpath(output_dir, "props_$(stub)_No_Disp_Lim.csv"), DataFrame(props, :auto))

            # Richness summary
            species_richness = length(unique(vec(comm)))
            row = (species_richness = species_richness,
                   sigma_str = sigma_str,
                   lambda_str = lambda_str,
                   neigh_str = neigh_str,
                   nsp_str = nsp_str,
                   case_label = cs.label)

            lock(result_lock) do
                push!(richness_results, row)
            end
        end)
    end

    foreach(fetch, tasks)
end

sleep(1)  # allow I/O to flush

# ============== Save final richness summary ==============
df_richness = DataFrame(richness_results)
CSV.write(joinpath(output_dir, "Richness_vs_nspecies.csv"), df_richness)


#using Pkg
#Pkg.activate(".")  # activates current directory as project
#cd("C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims")
#include("run/Figure_Regional_Pool_Size.jl")
