module HP_JCE_Sims

# Include submodules
include("Parameters.jl")
using .Parameters: SimulationConfig  

include("Helpers.jl")
include("SpeciesUtils.jl")
include("Environment.jl")
include("Dispersal.jl")
include("Competition.jl")
include("Simulation.jl")

# Export types and functions
export SimulationConfig, run_simulation

end
