using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# Choose which batch to run
include("batch_roughness.jl")  # Switch this line to other batch files as needed
