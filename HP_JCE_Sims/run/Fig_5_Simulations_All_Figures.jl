using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include("../src/HP_JCE_Sims.jl")
using .HP_JCE_Sims
using .HP_JCE_Sims.Simulation
using .HP_JCE_Sims.Parameters
using .HP_JCE_Sims.Environment: generate_habitat
using .HP_JCE_Sims.Simulation: get_neighbors, mean_kappa_M

using CSV, DataFrames
using Base.Threads: @spawn, @sync
using Random, Statistics

# === Output directory ===
base_output_dir = joinpath(@__DIR__, "..", "Outputs", "Fig_5")
isdir(base_output_dir) || mkpath(base_output_dir)

# === Shared simulation parameters ===
grid_size = 175
n_species = 500
n_generations = 50000
Moore_Neigh = 3

# === Threaded simulation runner ===
function run_simulations_set(sim_configs::Vector{<:NamedTuple}, subfolder::String)
    output_dir = joinpath(base_output_dir, subfolder)
    isdir(output_dir) || mkpath(output_dir)

    results = Vector{NamedTuple{(:species_richness, :param_str, :mean_kappa), Tuple{Int, String, Float64}}}()
    result_lock = Threads.SpinLock()
    tasks = []

    for cfg in sim_configs
        push!(tasks, @spawn begin
            config = SimulationConfig(
                grid_size = grid_size,
                n_species = n_species,
                n_generations = n_generations,
                sigma = cfg.sigma,
                JC_strength = cfg.JC_strength,
                dispersal_on = cfg.dispersal_on,
                dispersal_limit = 7,
                mortality_rate = 0.2,
                immigration_rate = .3*1e-5,
                dispersal_scale = cfg.dispersal_scale,
                Moore_Neigh = Moore_Neigh,
                spatial_autocorr = true,
                roughness_lambda = cfg.roughness_lambda,
                range_param_sim = 0.15,
                fecundity_meanlog = 0.0,
                fecundity_sdlog = 0.3,
                temporal_rho = 0.0
            )

            comm, props = run_simulation(config)

            mean_kappa = mean_kappa_M(comm, Moore_Neigh; threshold = 25)

            param_str = cfg.name

            CSV.write(joinpath(output_dir, "comm_$(param_str).csv"), DataFrame(comm, :auto))
            CSV.write(joinpath(output_dir, "props_$(param_str).csv"), DataFrame(props, :auto))

            richness = length(unique(vec(comm)))
            result = (species_richness = richness, param_str = param_str, mean_kappa = mean_kappa)

            lock(result_lock) do
                push!(results, result)
            end
        end)
    end

    @sync begin
        foreach(fetch, tasks)
    end

    sleep(1)

    df = DataFrame(results)
    CSV.write(joinpath(output_dir, "Summary_$(subfolder).csv"), df)
end

# === Simulation Sets ===

# Set 1: Vary JC_strength, fixed sigma = 99999, dispersal_scale = 1, dispersal_on = true
jc_strength_vals = exp10.(range(log10(0.1), log10(100), length=6))
sim_set_1 = [
    (; name = "JC_$(round(a, digits=3))", JC_strength = a, sigma = 99999.0, dispersal_on = true, dispersal_scale = 1.0, roughness_lambda = 0.1)
    for a in jc_strength_vals
]

# Set 2: Vary dispersal_scale, fixed JC_strength = 0.5, sigma = 99999
D_vals = [1.00, 2.75, 4.50, 6.25, 8.00, "none"]
sim_set_2 = [
    (; name = D == "none" ? "Disp_none" : "Disp_$(string(D))", JC_strength = 0.5, sigma = 99999.0,
       dispersal_on = D != "none", dispersal_scale = D == "none" ? 1.0 : D, roughness_lambda = 0.1)
    for D in D_vals
]

# Set 3: Vary roughness_lambda, fixed JC_strength = 0, sigma = 0.05, dispersal_off
lambda_vals =  exp10.(range(log10(1/300), log10(30.0), length=6))

sim_set_3 = [
    (; name = "Rough_$(replace(string(λ), "." => "_"))", JC_strength = 0.0, sigma = 0.05,
       dispersal_on = false, dispersal_scale = 1.0, roughness_lambda = λ)
    for λ in lambda_vals
]

# Set 4: Same as above, but JC_strength = 0.05
sim_set_4 = [
    (; name = "Rough_$(replace(string(λ), "." => "_"))", JC_strength = 0.5, sigma = 0.05,
       dispersal_on = false, dispersal_scale = 1.0, roughness_lambda = λ)
    for λ in lambda_vals
]

# === Run all sets ===
@sync begin
    tasks = [
        @spawn run_simulations_set(sim_set_1, "SimSet1_JCstrength"),
        @spawn run_simulations_set(sim_set_2, "SimSet2_Dispersal"),
        @spawn run_simulations_set(sim_set_3, "SimSet3_HP_only"),
        @spawn run_simulations_set(sim_set_4, "SimSet4_HP_plus_JCE")
    ]
    foreach(fetch, tasks)
end



# using Pkg
#Pkg.activate(".")  # activates current directory as project
#cd("C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims")
#include("run/Fig_5_Simulations_All_Figures.jl")