module BatchRuns

using ..Simulation
using ..Parameters

export run_multiple_simulations

"""
    run_multiple_simulations(configs::Vector{SimulationConfig})

Runs multiple simulations in parallel using multithreading.
Returns a vector of (community, proportions) tuples.
"""
function run_multiple_simulations(configs::Vector{SimulationConfig})
    n = length(configs)
    results = Vector{Tuple{Matrix{Int}, Matrix{Float64}}}(undef, n)

    Threads.@threads for i in 1:n
        results[i] = run_simulation(configs[i])
    end

    return results
end

end
