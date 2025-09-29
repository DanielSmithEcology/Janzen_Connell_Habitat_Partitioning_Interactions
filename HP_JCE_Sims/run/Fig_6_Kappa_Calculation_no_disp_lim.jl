#!/usr/bin/env julia
# ==========================================================
# Script: Fig_6_Kappa_Calculation.jl
# ==========================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include("../src/HP_JCE_Sims.jl")
using .HP_JCE_Sims
using .HP_JCE_Sims.Simulation
using .HP_JCE_Sims.Parameters
using .HP_JCE_Sims.Environment: generate_habitat
using .HP_JCE_Sims.Simulation: scaled_JCE_habitat_covariance

using CSV, DataFrames, Random
using Statistics                     # <-- mean, var
using Base.Threads: @spawn, @sync, nthreads

# ========= I/O =========
base_dir = joinpath(@__DIR__, "..", "Outputs", "Fig_6", "simulation_outputs_Fig_6", "Fig_6_JCE_Only")
isdir(base_dir) || error("Directory not found: ", base_dir)

out_path = joinpath(base_dir, "Summary_JCE_Only_Kappa_Values_no_disp_lim.csv")

# ========= Helpers =========

"Parse M from filename; GLOBAL -> 99"
function parse_M_from_filename(fname::AbstractString)::Int
    m = match(r"_M(\d+)_", fname)
    if m !== nothing
        return parse(Int, m.captures[1])
    elseif occursin("GLOBAL", fname)
        return 99
    else
        error("Could not parse M from filename: $fname")
    end
end

"Robustly read community matrix (drop non-numeric cols; round to Int)."
function read_community_matrix(fpath::AbstractString)::Matrix{Int}
    df = CSV.read(fpath, DataFrame)
    numcols = [c for c in names(df) if eltype(df[!, c]) <: Real]
    isempty(numcols) && error("No numeric columns in $fpath")
    A = Matrix{Float64}(df[:, numcols])
    A[.!isfinite.(A)] .= 0.0
    return round.(Int, A)
end

"Toroidal M×M Moore neighborhood counts via circshift."
function local_counts_toroidal(presence::BitMatrix, M::Int)::Matrix{Int}
    half = (M - 1) ÷ 2
    acc = zeros(Int, size(presence))
    @inbounds for di in -half:half, dj in -half:half
        acc .+= Int.(circshift(presence, (di, dj)))
    end
    return acc
end

"""
Average κ across species with abundance ≥ threshold.
κ_{i,M} = var(n_{i,M}) / mean(n_{i,M}), with toroidal neighborhoods.
"""
function mean_kappa_M(community::Matrix{Int}, M::Int; threshold::Int=10)::Float64
    sp = filter(!=(0), unique(vec(community)))
    kappas = Float64[]
    for s in sp
        presence = community .== s
        abundance = count(identity, presence)
        abundance < threshold && continue

        lc = local_counts_toroidal(presence, M)
        μ  = mean(vec(lc))
        μ == 0 && continue
        σ2 = var(vec(lc))
        push!(kappas, σ2 / μ)
    end
    return isempty(kappas) ? NaN : mean(kappas)
end

# ========= Discover files & choose window for GLOBAL =========
files = filter(fn -> endswith(lowercase(fn), ".csv") && startswith(basename(fn), "comm_"),
               readdir(base_dir; join=true))
isempty(files) && error("No comm_*.csv files found in: ", base_dir)

Ms = [parse_M_from_filename(f) for f in files]
finite_Ms = filter(!=(99), Ms)
isempty(finite_Ms) && error("No finite M files present; cannot set GLOBAL window.")
max_finite_M = maximum(finite_Ms)

M_calc_for(Mp::Int) = (Mp == 99 ? max_finite_M : Mp)

@info "Found $(length(files)) files. Threads: $(nthreads()). Largest finite M = $max_finite_M."

# ========= Compute (multi-threaded) =========
results = Vector{NamedTuple{(:M, :Kappa)}}(undef, length(files))

@sync for (idx, fpath) in enumerate(files)
    @spawn begin
        Mp    = parse_M_from_filename(fpath)
        Mcalc = M_calc_for(Mp)
        if iseven(Mcalc)
            @warn "Window M=$Mcalc from $fpath is even; incrementing to $(Mcalc+1)."
            Mcalc += 1
        end

        comm = read_community_matrix(fpath)

        # Optional: quick diagnostic of species passing threshold
        thr = 10
        # n_passing = sum(count(==(s), vec(comm)) ≥ thr for s in unique(vec(comm)) if s != 0)
        # @info "File $(basename(fpath)): species ≥ $thr -> $n_passing"

        kappa = try
            mean_kappa_M(comm, Mcalc; threshold=thr)
        catch e
            @warn "Failed on $fpath with error: $e"
            NaN
        end

        results[idx] = (M = Mp, Kappa = kappa)
        @info "Done: $(basename(fpath))  =>  M=$(Mp) (calc window $(Mcalc)),  κ=$(kappa)"
    end
end

# ========= Save =========
sorted = sort!(collect(results), by = x -> x.M)
df_out = DataFrame(sorted)
CSV.write(out_path, df_out)
println("Saved: ", out_path)
