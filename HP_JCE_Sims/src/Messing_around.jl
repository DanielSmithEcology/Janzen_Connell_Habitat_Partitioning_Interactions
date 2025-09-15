include("src/Dispersal.jl")  # or include your top-level module

include("src/HP_JCE_Sims.jl")        # reload the whole module

shape = 8.0
radius = 7
kernel, D_Lim = HP_JCE_Sims.Dispersal.dispersal_vector(shape, radius)

println("Kernel sum = ", sum(kernel))         # should be 1.0
println("Local dispersal proportion = ", D_Lim)

# Reshape for visualization
k_mat = reshape(kernel, (2*radius+1, 2*radius+1))


cd("C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims")
using Pkg
Pkg.activate(".")
include("src/HP_JCE_Sims.jl")

using .HP_JCE_Sims
using .HP_JCE_Sims.Parameters
using .HP_JCE_Sims.Simulation
using .HP_JCE_Sims.SpeciesUtils


config = SimulationConfig(
    grid_size = 125,
    n_species = 300,
    n_generations = 2500,
    sigma = 1.0,
    JC_strength = 1.0,
    dispersal_on = true,
    dispersal_limit = 6,
    mortality_rate = 0.2,
    immigration_rate = 0.00001,
    dispersal_scale = 5.0,
    Moore_Neigh = 1
)

comm, props = run_simulation(config)
@show size(comm)
@show props


@time comm, props = run_simulation(config)


using Plots


# Plot richness

using CSV, DataFrames, Glob

output_dir = joinpath(@__DIR__,  "Outputs", "Fig_3", "simulation_outputs_Fig_3G")
#cd("C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims/Outputs/Fig_3/simulation_outputs_Fig_3G")

output_dir = joinpath(@__DIR__, "Outputs", "Fig_3", "simulation_outputs_Fig_3G")
@show output_dir
@show isdir(output_dir)
@show readdir(output_dir)

using Glob

function species_richness_over_time(times::AbstractVector{Int}, p::AbstractMatrix{<:Real})
    richness = Int[]
    for t in times
        present = count(x -> x > 0, p[t, :])
        push!(richness, present)  
    end
    return richness
end


using Glob

output_dir = "C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims/Outputs/Fig_3/simulation_outputs_Fig_3G"


# go to project root (parent of `src`)
project_root = normpath(joinpath(@__DIR__, ".."))

output_dir = joinpath(project_root, "Outputs", "Fig_3", "simulation_outputs_Fig_3C")

out_norm = replace(output_dir, '\\' => '/')
props = Glob.glob(joinpath(out_norm, "props_*.csv"))


#props = Glob.glob(joinpath(output_dir, "props_*.csv"))
comm_files = Glob.glob(joinpath(output_dir, "props_*.csv"))


comm_files = filter(f ->
    startswith(basename(f), "props_") &&
    endswith(lowercase(f), ".csv"),
    readdir(output_dir; join=true)
)


#comm_files = Glob.glob(joinpath(output_dir, "props_*.csv"))
@show length(comm_files) comm_files

# Get all community files
comm_files = Glob.glob("props_*.csv", output_dir)

props = CSV.read(comm_files[21], DataFrame)
times = 1:size(props, 1)

comm_files[21]

println(comm_files)

using Plots
xx = 16
props = CSV.read(comm_files[xx], DataFrame)
comm_files[xx]

richness = species_richness_over_time(times, Matrix(props))
plot((times), richness, xlabel="Time", ylabel="Richness", label="", lw=2, title="Species Richness Over Time")
#plot(times, richness, xlabel="Time", ylabel="Richness", label="", lw=2, title="Species Richness Over Time")

# Plot species proportions
plot_species_proportions(times, props)


using Profile
@profile comm, props = run_simulation(config)
Profile.print()



cd("C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims")
using Pkg
Pkg.activate(".")

using Random  # ← Add this
Random.seed!(1234)  # ← Set your reproducible seed

include("src/HP_JCE_Sims.jl")

using .HP_JCE_Sims
using .HP_JCE_Sims.Parameters
using .HP_JCE_Sims.Simulation
using .HP_JCE_Sims.SpeciesUtils
using .HP_JCE_Sims.Environment
using .HP_JCE_Sims.Competition
using .HP_JCE_Sims.Dispersal

# Define configuration with spatial autocorrelation
config = SimulationConfig(
    grid_size = 125,
    n_species = 300,
    n_generations = 5000,
    sigma = .025,
    JC_strength = 1.0,
    dispersal_on = false,
    dispersal_limit = 6,
    mortality_rate = 0.2,
    immigration_rate = 0.00001,
    dispersal_scale = 5.0,
    Moore_Neigh = 2,
    spatial_autocorr = true,
    kernel_sigma = 5.0,
    roughness_lambda = 0.0  # smooth surface with some noise
)

#comm, props = run_simulation(config)

@time comm, props = run_simulation(config)



@show size(comm)
@show props


using Plots
times = 1:size(Props, 1)

function species_richness_over_time(times::AbstractVector{Int}, p::AbstractMatrix{<:Real})
    richness = Int[]
    for t in times
        present = count(x -> x > 0, p[t, :])
        push!(richness, present)  
    end
    return richness
end


# Plot richness
richness = species_richness_over_time(times, Matrix(Props))
plot(times, richness, xlabel="Time", ylabel="Richness", label="", lw=2, title="Species Richness Over Time")






using ImageFiltering
using Random
using Statistics
using LinearAlgebra

# ========== Function to test ==========
"""
    update_dispersal_pressure!(
        dispersal_pressure, community, dispersal_kernel_2d, p, D_Lim
    )

Compute local dispersal for each species and add global fallback.
"""
function update_dispersal_pressurex!(
    dispersal_pressure::Array{Float64, 3},
    community::Matrix{Int},
    dispersal_kernel_2d::Matrix{Float64},
    p::Vector{Float64},
    D_Lim::Float64
)
    present_species = filter(!=(0), unique(vec(community)))

    for s in present_species
        presence_map = @. community == s
        local_disp = imfilter(Float64.(presence_map), dispersal_kernel_2d, Pad(:circular))
        dispersal_pressure[:, :, s] .= D_Lim * local_disp .+ (1 - D_Lim) * p[s]
    end

    # Optional: handle extinct species
    extinct_species = setdiff(1:length(p), present_species)
    for s in extinct_species
        dispersal_pressure[:, :, s] .= (1 - D_Lim) * p[s]
    end
end

# ========== Test Setup ==========
G = 25         # grid size
S = 3          # number of species
r = 6          # dispersal radius
D_Lim = 0.8    # local dispersal strength

# Make a dummy dispersal kernel (Chebyshev shape)
function make_kernel(r)
    k = zeros(2r+1, 2r+1)
    for dx in -r:r, dy in -r:r
        d = max(abs(dx), abs(dy))
        k[dx+r+1, dy+r+1] = 1 / (1 + d)
    end
    return k ./ sum(k)
end

disp_kernel_2d = make_kernel(r)

# Dummy community with 3 species (1, 2, 3) and some empty cells (0)
community = fill(0, G, G)
community[10:12, 10:12] .= 1
community[5:6, 18:19] .= 2
community[20:21, 5:6] .= 3

# Regional proportions (can be anything summing to 1)
p = [0.3, 0.4, 0.3]

# Output array
dispersal_pressure = zeros(G, G, S)

# ========== Run Test ==========
update_dispersal_pressurex!(dispersal_pressure, community, disp_kernel_2d, p, D_Lim)

# ========== Inspect Output ==========
using CairoMakie

fig = Figure(size=(800, 250))
for s in 1:S
    ax = Axis(fig[1, s], title="Species $s")
    CairoMakie.heatmap!(ax, dispersal_pressure[:, :, 1])
end
fig

Plot.heatmap(dispersal_pressure[:, :, 1])

using Plots
Plots.heatmap(dispersal_pressure[:, :, 3])




using CairoMakie
using Statistics

# === Dispersal kernel logic (from your real simulation code) ===
"""
    dispersal_vector(scale::Float64, r::Int) -> (Vector, Float64)

Returns the 1D dispersal vector (Chebyshev-based) and D_Lim.
"""
function dispersal_vectorx(scale::Float64, r::Int)
    distances = 0:r
    weights = exp.(-distances ./ scale)
    return weights ./ sum(weights), sum(weights)
end


function dispersal_vectorx2(Az::Float64, maxlen::Int)
    dx = 1e-4
    disp_vec = Float64[]

    # First bin: from 0 to 0.5
    x0 = 0:dx:0.5
    y0 = dispersal_kernel.(x0, Az)
    bin_prob = integrate(x0, y0)
    push!(disp_vec, bin_prob)
    total_prob = bin_prob

    d = 1
    while d < maxlen && total_prob < 0.99
        x = (d - 0.5):dx:(d + 0.5)
        y = dispersal_kernel.(x, Az)
        bin_prob = (1 / (8 * d)) * integrate(x, y)

        push!(disp_vec, bin_prob)
        total_prob += bin_prob*(d*8)
        d += 1
    end

    return disp_vec, total_prob
end


"""
    create_2d_kernel(kernel::Vector) -> Matrix

Expands a 1D kernel into a 2D kernel assuming Chebyshev distance.
"""
function create_2d_kernelx(kernel::Vector{Float64})
    r = length(kernel) - 1
    K = zeros(Float64, 2r + 1, 2r + 1)
    for dx in -r:r, dy in -r:r
        d = max(abs(dx), abs(dy))
        K[dx + r + 1, dy + r + 1] = kernel[d + 1]
    end
    return K ./ sum(K)  # normalize
end

"""
    update_dispersal_pressurex!(
        dispersal_pressure, community, dispersal_kernel_2d, p, D_Lim
    )
"""
function update_dispersal_pressurex!(
    dispersal_pressure::Array{Float64, 3},
    community::Matrix{Int},
    dispersal_kernel_2d::Matrix{Float64},
    p::Vector{Float64},
    D_Lim::Float64
)
    G = size(community, 1)
    S = length(p)

    present_species = filter(!=(0), unique(vec(community)))

    for s in present_species
        presence_map = @. community == s
        local_disp = imfilter(Float64.(presence_map), dispersal_kernel_2d, Pad(:circular))
        dispersal_pressure[:, :, s] .= D_Lim * local_disp .+ (1 - D_Lim) * p[s]
    end

    # Handle extinct species
    extinct_species = setdiff(1:S, present_species)
    for s in extinct_species
        dispersal_pressure[:, :, s] .= (1 - D_Lim) * p[s]
    end
end


function dispersal_kernel(x::Real, Az::Real)
    return (2 * π * x / (π * Az^2)) * exp(-x^2 / Az^2)
end

# === Parameters ===
G = 25
S = 3
r = 7
scale = 1.5

using Plots
# === Create kernel ===

using NumericalIntegration


disp_vector_1d, _ = dispersal_vectorx2(scale, r)
aaa = vec([ disp_vector_1d[1], disp_vector_1d[2]*8, disp_vector_1d[3]*16])

aaa = vec([ disp_vector_1d[1], disp_vector_1d[2]*8, disp_vector_1d[3]*16, disp_vector_1d[4]*24 ])
plot(aaa)
Plots.plot(disp_vector_1d)

disp_kernel_2d = create_2d_kernelx(disp_vector_1d)

# === Dummy community ===
community = fill(0, G, G)
community[10:12, 10:12] .= 1
community[5:6, 18:19] .= 2
community[20:21, 5:6] .= 3

# === Regional proportions ===
p = [0.3, 0.4, 0.3]

# === Output array ===
dispersal_pressure = zeros(G, G, S)

# === Run update ===
update_dispersal_pressurex!(dispersal_pressure, community, disp_kernel_2d, p, 0.8)

# === Plot results ===
fig = Figure(size=(900, 300))
for s in 1:S
    ax = Axis(fig[1, s], title="Dispersal Pressure: Species $s")
    heatmap!(ax, dispersal_pressure[:, :, s])
end
fig
using Plots
Plots.heatmap(dispersal_pressure[:, :, 3])



using ImageFiltering
function generate_habitat_with_noisex(G::Int; σ::Float64=5.0, λ::Float64=0.1)
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

G = 100
σ = 5.0
λ = .1
fff = generate_habitat_with_noisex(G; σ=σ, λ=λ)

using Plots
heatmap(fff)

histogram(vec(fff), bins = 30, color = :skyblue)



using ImageFiltering, Plots

function generate_habitat_with_noisexx(G::Int; σ::Float64=5.0, λ::Float64=0.1)
    base = rand(G, G)
    kernel = Kernel.gaussian((σ, σ))
    smoothed = imfilter(base, kernel, Pad(:circular))
    blended = λ * base .+ (1 - λ) * smoothed

    # Match distribution back to original base
    sorted_blended = sort(vec(blended))
    sorted_base = sort(vec(base))
    flat_blended = vec(blended)
    ranks = [searchsortedfirst(sorted_blended, x) for x in flat_blended]
    matched = reshape(sorted_base[ranks], G, G)

    return matched
end

using Random, Statistics
Random.seed!(42)
# Test it
G = 175
σ = 10.0
λ = 0.0
fff = generate_habitat_with_noisexx(G; σ=σ, λ=λ)
heatmap(fff,aspect_ratio=1)
histogram(vec(fff))
fff0 = copy(fff)

compute_morans_I(fff)

using SpatialDependence
W = contiguity_weights(G, "queen")  # define spatial weights on grid cells
I = morans_I(vec(habitat_matrix), W)



using Pkg
Pkg.add("SpatialDependence")

using Pkg
Pkg.add("GeoStatsFunctions")
Pkg.add("GeoStats")

using GeoStatsFunctions
using GeoStats

# Suppose habitat is a G×G matrix
X, Y = collect(1:G), collect(1:G)
coords = [(i, j) for i in X, j in Y] |> vec
values = vec(habitat)

# Build point data
df = georef((value=values,), coords)

ev = EmpiricalVariogram(df, :value; maxlag=G/2, nlags=20)



function compute_morans_I(matrix::AbstractMatrix)
    G = size(matrix, 1)
    N = G * G
    x̄ = mean(matrix)
    diffs = matrix .- x̄
    s2 = sum(diffs.^2)

    # Moore neighborhood (excluding self), using mod1 for wrap-around
    function neighbors(i, j)
        [(mod1(i + di, G), mod1(j + dj, G))
         for di in -1:1, dj in -1:1 if !(di == 0 && dj == 0)]
    end

    w_sum = 0.0
    num = 0.0

    for i in 1:G, j in 1:G
        xi = matrix[i, j] - x̄
        for (ni, nj) in neighbors(i, j)
            xj = matrix[ni, nj] - x̄
            num += xi * xj
            w_sum += 1
        end
    end

    return (N / w_sum) * (num / s2)
end


Pkg.add("GaussianRandomFields")

using GaussianRandomFields, Plots, GeoStatsFunctions
ρ = 1/30
ν = 5/4
pts = range(0, stop=1, step=1/174)

ffs = rand(175,175)

cov = CovarianceFunction(2, Matern(ρ, ν))
grf = GaussianRandomField(cov, CirculantEmbedding(),  pts, pts, minpadding=625)
field = sample(grf)
heatmap(field)
histogram(vec(field))

f2 = match_distribution(field,ffs)
histogram(vec(f2))
heatmap(f2)
# 2. Empirical variogram
ev = EmpiricalVariogram(field)



for ν in [5/4, 1, 7/8, 3/4, 5/8, 1/2]
    cov = CovarianceFunction(2, Matern(1/4, ν))
    grf = GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=625)
    heatmap(grf)
end





using Random, Statistics
using Plots
using ImageFiltering

# === Dummy versions of your existing functions ===

function gaussian_smooth(mat::Matrix{Float64}, sigma::Float64)
    kernel = Kernel.gaussian(sigma)
    return imfilter(mat, kernel; border="replicate")
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

function generate_raw_noise(rows::Int, cols::Int, lambda::Float64)
    return randn(rows, cols) .* lambda
end

function update_habitat!(
    current_hab::Matrix{Float64},
    target_dist::Matrix{Float64},
    kernel_sigma::Float64,
    lambda::Float64,
    rho::Float64
)
    new_noise = generate_raw_noise(size(current_hab)..., lambda)
    smoothed_new = gaussian_smooth(new_noise, kernel_sigma)
    matched_new = match_distribution(smoothed_new, target_dist)
    current_hab .= rho .* current_hab .+ (1.0 - rho) .* matched_new
end

# === Simulation parameters ===

G = 175
T = 10
lambda = 1.0
kernel_sigma = 1.5
rho = 0.8  # temporal autocorrelation

# === Initial habitat ===

heatmap(fff,aspect_ratio=1)

Random.seed!(42)
init_noise = generate_raw_noise(G, G, lambda)
habitat0 = gaussian_smooth(init_noise, kernel_sigma)
habitat0 = match_distribution(habitat0, habitat0)  # Normalize
habitat = copy(habitat0)

# === Plot over time ===

plot(layout = (2, 5), size=(1200, 500), title="Temporal Habitat Dynamics", dpi=150)

for t in 1:T
    update_habitat!(habitat, habitat0, kernel_sigma, lambda, rho)
    heatmap!(habitat, title="t = $t", clims=(minimum(habitat0), maximum(habitat0)))
end

display(current())




using ImageFiltering

# 1. Gaussian smoothing
function gaussian_smooth(mat::Matrix{Float64}, sigma::Float64)
    kernel = Kernel.gaussian((sigma, sigma))
    return imfilter(mat, kernel; border="replicate")
end

using ImageFiltering

using ImageFiltering

function gaussian_smooth(mat::Matrix{Float64}, sigma::Float64)
    kernel = Kernel.gaussian((sigma, sigma))
    return imfilter(mat, kernel, Pad(:replicate))  # use Pad as a positional argument
end


# 2. Match distribution of source to target (preserve marginal)
function match_distribution(source::Matrix{<:Real}, target::Matrix{<:Real})
    flat_source = vec(source)
    flat_target = vec(target)

    sorted_source = sort(flat_source)
    sorted_target = sort(flat_target)

    # Compute rank order and remap
    ranks = sortperm(flat_source)
    matched = similar(flat_source)
    matched[ranks] = sorted_target

    return reshape(matched, size(source))
end

# 3. Additive noise to simulate habitat fluctuations
function generate_raw_noise(rows::Int, cols::Int, lambda::Float64)
    return randn(rows, cols) .* lambda
end

# 4. Temporal update of habitat with controlled smooth noise
function update_habitat_temporally!(
    habitat::Matrix{Float64},
    kernel_sigma::Float64,
    noise_lambda::Float64
)
    noise = generate_raw_noise(size(habitat)..., noise_lambda)
    smoothed_noise = gaussian_smooth(noise, kernel_sigma)
    updated = habitat .+ smoothed_noise
    matched = match_distribution(updated, habitat)
    habitat .= matched  # in-place update
end




function update_habitatx!(
    current_hab::Matrix{Float64},
    target_dist::Matrix{Float64},
    kernel_sigma::Float64,
    lambda::Float64,
    rho::Float64
)
    # New noise field
    new_field = generate_raw_noise(size(current_hab)..., lambda)

    # Smooth it
    smoothed_new = gaussian_smooth(new_field, kernel_sigma)

    # Match distribution of original target
    matched_new = match_distribution(smoothed_new, target_dist)

    # Blend with previous environment
    updated = rho .* current_hab .+ (1.0 - rho) .* matched_new

    current_hab .= updated
end


G = 100
initial_habitat = rand(G, G)
kernel_sigma = 2.0
noise_lambda = 0.05

for t in 1:10
    update_habitat_temporally!(initial_habitat, kernel_sigma, noise_lambda)
    println("Time $t updated.")
end

heatmap(fff)

G = 50
T = 10
lambda = 1.0
kernel_sigma = 1.5
rho = 0.9  # temporal autocorrelation

G = 100
kernel_sigma = 5.0
lambda = 0.05

Random.seed!(42)
# Test it
G = 175
σ = 10.0
λ = 0.0
fff = generate_habitat_with_noisexx(G; σ=σ, λ=λ)
heatmap(fff)

fff0 = copy(fff)



update_habitatx!(fff, fff0,kernel_sigma, lambda,rho)

heatmap(fff)
heatmap(fff0)
histogram(vec(fff))
histogram(vec(fff0))

initial_habitat


using NeutralLandscapes

siz = 200, 200



function match_distribution(source::Matrix{<:Real}, target::Matrix{<:Real})
    flat_source = vec(source)
    flat_target = vec(target)

    sorted_source = sort(flat_source)
    sorted_target = sort(flat_target)

    # Compute rank order and remap
    ranks = sortperm(flat_source)
    matched = similar(flat_source)
    matched[ranks] = sorted_target

    return reshape(matched, size(source))
end


baseline = rand(175,175)


Baseline = rand(MidpointDisplacement(0.1), siz)
histogram(vec(Baseline))

Fig1f = rand(MidpointDisplacement(0.99), siz)
Updated_Hab = match_distribution(Fig1f,baseline)
heatmap(Updated_Hab)
histogram(vec(Updated_Hab))


SpatiallyAutocorrelatedUpdater!(spatialupdater=MidpointDisplacement(0.5))


heatmap(Fig1f)
histogram(vec(Fig1f))

Updated_Hab = match_distribution(Fig1f,baseline)

heatmap(Updated_Hab)
histogram(vec(Updated_Hab))



heatmap(Fig1f)
histogram(vec(Fig1f))




Baseline = rand(MidpointDisplacement(0.1), siz)

using NeutralLandscapes

siz = 175,175

Random.seed!(154)
baseline = rand(175,175)
#histogram(vec(baseline))
Fig1f = rand(MidpointDisplacement(0.9999999999999), siz)
Fig1fb = match_distribution(Fig1f,baseline)
heatmap(Fig1fb)

histogram(vec(Fig1fb))


siz = 200, 200
baseline = rand(200,200)
Fig1f = rand(MidpointDisplacement(0.01), siz)
update!(SpatiallyAutocorrelatedUpdater(MidpointDisplacement(.999999999999),0,.01), Updated_Hab)
Updated_Hab = match_distribution(Updated_Hab,baseline)



heatmap(Fig1f)
Updated_Hab = match_distribution(Fig1f,baseline)
heatmap(Updated_Hab)

update!(SpatiallyAutocorrelatedUpdater(MidpointDisplacement(.999999999999),0,.01), Updated_Hab)
Updated_Hab = match_distribution(Updated_Hab,baseline)

heatmap(Updated_Hab)
histogram(vec(Updated_Hab))



"""
    compute_morans_I(matrix::AbstractMatrix)

Compute Moran's I for a 2D spatial matrix. Assumes toroidal (wrap-around) adjacency.
"""

"""
    compute_morans_I(matrix::AbstractMatrix)

Compute Moran's I for a 2D spatial matrix using 8-neighbor Moore neighborhood and toroidal wrapping.
"""
function compute_morans_I(matrix::AbstractMatrix)
    G = size(matrix, 1)
    N = G * G
    x̄ = mean(matrix)
    diffs = matrix .- x̄
    s2 = sum(diffs.^2)

    # Moore neighborhood (excluding self), using mod1 for wrap-around
    function neighbors(i, j)
        [(mod1(i + di, G), mod1(j + dj, G))
         for di in -1:1, dj in -1:1 if !(di == 0 && dj == 0)]
    end

    w_sum = 0.0
    num = 0.0

    for i in 1:G, j in 1:G
        xi = matrix[i, j] - x̄
        for (ni, nj) in neighbors(i, j)
            xj = matrix[ni, nj] - x̄
            num += xi * xj
            w_sum += 1
        end
    end

    return (N / w_sum) * (num / s2)
end




compute_morans_I(Updated_Hab)


function generate_habitat_with_noisexx(G::Int; σ::Float64=5.0, λ::Float64=0.1)
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

# Test it
Random.seed!(1234)
G = 200
σ = 8.0
λ = 0.04
fff = generate_habitat_with_noisexx(G; σ=σ, λ=λ)
heatmap(fff)
histogram(vec(fff))
compute_morans_I(fff)


compute_morans_I(fff)





sau = SpatiallyAutocorrelatedUpdater(spatialupdater=MidpointDisplacement(0.5))
makeplot(sau)


function makeplot(updater)
    init = rand(PlanarGradient(π/2), 100, 100)
    seq = update(updater, init, 30)
    print(length(seq))
    plot(heatmap(init, cbar=:none), heatmap(seq[10], cbar=:none), heatmap(seq[20], cbar=:none), heatmap(seq[30], cbar=:none), size=(700,200), layout=(1,4))
end

sau = SpatiallyAutocorrelatedUpdater(spatialupdater=MidpointDisplacement(0.5))
makeplot(sau)




stau = SpatiotemporallyAutocorrelatedUpdater(spatialupdater=MidpointDisplacement(0.5), variability=1.0)
makeplot(stau)
stau








using CSV, DataFrames
using Glob  # if not installed: Pkg.add("Glob")

output_dir = joinpath(@__DIR__,  "Outputs", "Fig_3", "simulation_outputs_Fig_3C")


# Get all community files
comm_files = Glob.glob("comm_*.csv", output_dir)

results = NamedTuple[]

for file in comm_files
    df = CSV.read(file, DataFrame)
    richness = length(unique(vec(Matrix(df))))

    # Extract σ, λ, and neigh from filename
    m = match(r"comm_neigh_(.+?)_sigma_(.+?)_lambda_(.+?)_No_Disp_Lim", basename(file))
    if m !== nothing
        neigh_str, σ_str, λ_str = m.captures
        push!(results, (
            species_richness = richness,
            σ_str = σ_str,
            λ_str = λ_str,
            neigh_str = neigh_str
        ))
    end
end

df_richness = DataFrame(results)
CSV.write(joinpath(output_dir, "Richness_M_sigmah_No_Disp_Lim_High_Autocorr_2.csv"), df_richness)


output_dir = joinpath(@__DIR__,  "Outputs", "Fig_3", "simulation_outputs_Fig_3C")
all_files = readdir(output_dir)

println("Files in output dir:")
foreach(println, all_files)



fs = CSV.read(comm_files[1],DataFrame)
unique(vec(Matrix(fs)))


function get_neighbors(community::Matrix{Int}, i::Int, j::Int, r::Int)
    G = size(community, 1)
    neighbors = Int[]
    for dx in -r:r, dy in -r:r
        # Uncomment this if you want to exclude the focal cell:
        # if dx == 0 && dy == 0
        #     continue
        # end
        ni = mod1(i + dx, G)
        nj = mod1(j + dy, G)
        push!(neighbors, community[ni, nj])
    end
    return neighbors
end
community = reshape(1:25, 5, 5)  # 5x5 grid, values 1 to 25

i, j = 3, 3  # center cell in 5x5 grid
for r in 0:2
    neighbors = get_neighbors(community, i, j, r)
    println("r = $r: ", neighbors)
end



function get_neighbors(community::Matrix{Int}, i::Int, j::Int, r::Int)
    G = size(community, 1)
    neighbors = Int[]
    for dx in -r:r, dy in -r:r
        ni = mod1(i + dx, G)
        nj = mod1(j + dy, G)
        push!(neighbors, community[ni, nj])
    end
    return neighbors
end

community = collect(reshape(1:25, 5, 5))

i, j = 3, 3  # center cell
for r in 0:2
    println("r = $r: ", get_neighbors(community, i, j, r))
end




neigh =0

JanCon_a, neigh_val = neigh == "none" ? (0.0, 0) : (0.5, neigh) # if no JCEs, set a =0.0; otherwise set to 0.5


σ = 0.015
replace(string(round(σ, digits=4)), "." => "_")



using GaussianRandomFields

# Define your grid
G = 100
x = range(0, stop=1, length=G)
y = range(0, stop=1, length=G)

# Define covariance structure (approximating your equation)
r = 0.2  # range
s = 1.0  # sill
n = 0.05 # nugget

covfun = CovarianceFunction(2, Gaussian(s))  # This uses exp(- (d/r)^2)
grf = GaussianRandomField(covfun, CirculantEmbedding(), x, y)
field = sample(grf)
heatmap(field)

using GeoStatsFunctions, Plots

# Compute empirical variogram (2D isotropic assumed)
empvario = EmpiricalVariogram(field)

# Plot
plot(empvario)


Pkg.add("Tables")
Pkg.add("GeoTables")


NEIGH = [1 2 3 4 4 5 6 7 8 8 8 9 9]
a=.2
S=10
optima         = rand(S)
σ = .05
hab=.3
habitat_effects = @. exp(-((hab - optima)^2) / (2 * σ^2))


@. exp(-a* count(==(s), NEIGH)) where s=1:S


jc_penalties = zeros(S)
for s in 1:S
    jc_penalties[s] = exp(-a * count(==(s), NEIGH))
end


using GaussianRandomFields
using GeoStatsFunctions
using Tables, GeoTables  # ✅ Required for GeoTable
using Plots

# === Generate GRF on grid ===
Random.seed!(154)

G = 175
sill = 1.0
range2 = 0.2
nugget = 0.00

scale = sqrt(3) / range2
x = range(0, 1, length=G) .* scale
y = range(0, 1, length=G) .* scale

covfun = CovarianceFunction(2, Gaussian(sill, range2))

grf = GaussianRandomField(covfun, CirculantEmbedding(), x, y)
field = sample(grf)
field_nugget = field .+ sqrt(nugget) .* randn(size(field))

heatmap(field_nugget)
baseline = rand(175,175)
Fig1fb = match_distribution(field_nugget,baseline)
heatmap(Fig1fb)


#histogram(vec(baseline))
Fig1f = rand(MidpointDisplacement(0.9999999999999), siz)
Fig1fb = match_distribution(Fig1f,baseline)
heatmap(Fig1fb)

histogram(vec(Fig1fb))


cov = CovarianceFunction(2, Matern(1/4, 3/4)) # length scale 1/4, smoothness 3/4

pts = range(0, stop=1, step=1/400)
grf = GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=113)
grf2 = match_distribution(grf,baseline)

heatmap(grf)
heatmap(grf2)






using GaussianRandomFields

function match_distribution(source::Matrix{<:Real}, target::Matrix{<:Real})
    flat_source = vec(source)
    flat_target = vec(target)

    sorted_source = sort(flat_source)
    sorted_target = sort(flat_target)

    # Compute rank order and remap
    ranks = sortperm(flat_source)
    matched = similar(flat_source)
    matched[ranks] = sorted_target

    return reshape(matched, size(source))
end


Random.seed!(154)
baseline = rand(175,175)

G = 175
sill = 1.0
range2 = 0.2         # this is 'r' in the paper
nugget = 0.0

# Adjust range so that exp(-3(d/r)^2) = exp(-(d/λ)^2)
λ = range2 / sqrt(3)

# Grid
x = range(0, 1, length=G)
y = range(0, 1, length=G)

# Covariance with adjusted range
covfun = CovarianceFunction(2, Gaussian(λ, σ = sqrt(sill)))
grf = GaussianRandomField(covfun, CirculantEmbedding(), x, y)
field = sample(grf)

# Add nugget effect
field_nugget = field .+ sqrt(nugget) .* randn(size(field))

Final_Habitat = match_distribution(field_nugget,baseline)



function rescale_to_uniform(field::Matrix{<:Real})
    flat = vec(field)
    ranks = sortperm(flat)
    matched = similar(flat)
    matched[ranks] = range(0, 1, length=length(flat))
    return reshape(matched, size(field))
end

Final_Habitat = rescale_to_uniform(field_nugget)

heatmap(Final_Habitat)




using GaussianRandomFields, Random

"""
    generate_habitat(G; sill=1.0, range=0.2, nugget=0.0, seed=nothing)

Generates a `G x G` habitat matrix as a Gaussian random field with optional nugget effect
and rescales the field to match a uniform distribution.

# Arguments
- `G`: grid size (e.g., 175)
- `sill`: controls the large-scale variance (default 1.0)
- `range`: spatial range parameter `r` (as in the paper)
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

generaG = 175
habitat = generate_habitat(G; sill=1.0, range_param=0.1, nugget=0.0, seed=42)
heatmap(habitat)







# === Flatten to coordinate-value pairs ===
coords = [(Float64(i), Float64(j)) for i in 1:G, j in 1:G]
values = vec(field_nugget)

# === Wrap as GeoTable ===
data = GeoTable((location=coords, stress=values))

# === Compute empirical variogram ===
empvario = EmpiricalVariogram(data, :stress)

# === Fit model ===
fit_model = fit(GaussianVariogram, empvario)

# === Plot ===
plot(empvario, label="Empirical", xlabel="Distance", ylabel="γ(d)")
plot!(fit_model, label="Gaussian Fit", lw=2)

# === Output estimated params ===
println("Estimated range: ", fit_model.range)
println("Estimated sill: ", fit_model.sill)
println("Estimated nugget: ", fit_model.nugget)



rr = rand(G,G)
rr2 = match_distribution(field_nugget,rr)
heatmap(rr2)
heatmap(field_nugget)











Random.seed!(154)
baseline = rand(175,175)
#histogram(vec(baseline))
Fig1f = rand(MidpointDisplacement(0.9999999999999), siz)
Fig1fb = match_distribution(Fig1f,baseline)
heatmap(Fig1fb)

histogram(vec(Fig1fb))


cov = CovarianceFunction(2, Matern(1/4, 3/4)) # length scale 1/4, smoothness 3/4

pts = range(0, stop=1, step=1/400)
grf = GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=113)
grf2 = match_distribution(grf,baseline)

heatmap(grf)
heatmap(grf2)







using GaussianRandomFields, Random

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

using GaussianRandomFields, Random

G=175
Habb = generate_habitat(G; sill=1.0, range_param=.15, nugget=0.0, seed=11)
heatmap(Habb)
Habb1 = generate_habitat(G; sill=1.0, range_param=.15, nugget=0.0, seed=11)
Habb2 = generate_habitat(G; sill=1.0, range_param=.15, nugget=1.0, seed=11)
Habb3 = generate_habitat(G; sill=1.0, range_param=.15, nugget=5.0, seed=11)
Habb4 = generate_habitat(G; sill=1.0, range_param=.15, nugget=0.0, seed=11)
Habb5 = generate_habitat(G; sill=1.0, range_param=.15, nugget=1.0, seed=11)
Habb6 = generate_habitat(G; sill=1.0, range_param=.15, nugget=5.0, seed=11)
Habb7 = generate_habitat(G; sill=1.0, range_param=.15, nugget=0.0, seed=11)
Habb8 = generate_habitat(G; sill=1.0, range_param=.15, nugget=1.0, seed=11)
Habb9 = generate_habitat(G; sill=1.0, range_param=.15, nugget=5.0, seed=11)


CSV.write


# Define output path (escaped backslashes or use raw string)
output_path = raw"C:\Users\smith\OneDrive\Desktop\CNDD New Ideas\Hab_Par_New\Code_For_Figures"

# Generate habitats
Habb1 = generate_habitat(G; sill=1.0, range_param=0.05, nugget=0.00, seed=11)
Habb2 = generate_habitat(G; sill=1.0, range_param=0.05, nugget=1.0, seed=11)
Habb3 = generate_habitat(G; sill=1.0, range_param=0.05, nugget=10.0, seed=11)
Habb4 = generate_habitat(G; sill=1.0, range_param=0.1, nugget=0.00, seed=11)
Habb5 = generate_habitat(G; sill=1.0, range_param=0.2, nugget=1.0, seed=11)
Habb6 = generate_habitat(G; sill=1.0, range_param=0.2, nugget=10.0, seed=11)
Habb7 = generate_habitat(G; sill=1.0, range_param=0.5,  nugget=0.00, seed=11)
Habb8 = generate_habitat(G; sill=1.0, range_param=0.5, nugget=1.0, seed=11)
Habb9 = generate_habitat(G; sill=1.0, range_param=0.5, nugget=10.0, seed=11)

# Convert to DataFrames and save
CSV.write(joinpath(output_path, "Habb1.csv"), DataFrame(Habb1, :auto); writeheader = false)
CSV.write(joinpath(output_path, "Habb2.csv"), DataFrame(Habb2, :auto); writeheader = false)
CSV.write(joinpath(output_path, "Habb3.csv"), DataFrame(Habb3, :auto); writeheader = false)
CSV.write(joinpath(output_path, "Habb4.csv"), DataFrame(Habb4, :auto); writeheader = false)
CSV.write(joinpath(output_path, "Habb5.csv"), DataFrame(Habb5, :auto); writeheader = false)
CSV.write(joinpath(output_path, "Habb6.csv"), DataFrame(Habb6, :auto); writeheader = false)
CSV.write(joinpath(output_path, "Habb7.csv"), DataFrame(Habb7, :auto); writeheader = false)
CSV.write(joinpath(output_path, "Habb8.csv"), DataFrame(Habb8, :auto); writeheader = false)
CSV.write(joinpath(output_path, "Habb9.csv"), DataFrame(Habb9, :auto); writeheader = false)




using Plots
heatmap(Habb9)
histogram(vec(Habb1))

using Statistics
compute_morans_I(Habb)

0.025, 0.05, 0.1, 0.15, 0.2, 0.5

0.0,  0.5, 1.0 , 2.5, 5.0 ,10.0

function get_neighborsx(community::Matrix{Int}, i::Int, j::Int, r::Int)
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

function mean_kappa_Mx(community::Matrix{Int}, M::Int; threshold::Int = 10)
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


NUM = 50         # maximum integer
rows = 175      # number of rows
cols = 175
# number of columns

mat = rand(1:NUM, rows, cols)


using Statistics
mean_kappa_Mx(mat, 3; threshold=25)
