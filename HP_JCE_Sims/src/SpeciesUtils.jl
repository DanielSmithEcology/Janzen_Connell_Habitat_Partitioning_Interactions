module SpeciesUtils
using Plots

export species_proportions, species_richness, species_richness_over_time, plot_species_proportions

"""
    species_proportions(n_species::Int, community::Matrix{Int})

Compute the proportion of each species (from 1 to `n_species`) in the community matrix.
Returns a vector of proportions.
"""
function species_proportions(n_species::Int, community::Matrix{Int})
    props = zeros(n_species)
    total = length(community)
    for i in 1:n_species
        props[i] = count(x -> x == i, community) / total
    end
    return props
end

"""
    species_richness(simulations::Vector)

Compute species richness from a list of simulations.
Each simulation is assumed to return a tuple (community, timeseries).
Returns a vector of richness values.
"""
function species_richness(simulations::Vector)
    richness = Int[]
    for sim in simulations
        community = sim[1]
        push!(richness, length(unique(community)))
    end
    return richness
end

"""
    species_richness_over_time(times::AbstractVector{Int}, p::AbstractMatrix{<:Real})

Compute the number of species present (nonzero) at each time step given a timeseries matrix `p`.
Each row of `p` is assumed to be a vector of species proportions.
"""
function species_richness_over_time(times::AbstractVector{Int}, p::AbstractMatrix{<:Real})
    richness = Int[]
    for t in times
        present = count(x -> x > 0, p[t, :])
        push!(richness, present)
    end
    return richness
end

function plot_species_proportions(times::AbstractVector{Int}, p::AbstractMatrix{<:Real})
    plot(times, p[:, 1], label = "Species 1", lw=2)
    for s in 2:size(p, 2)
        plot!(times, p[:, s], label = "Species $s", lw=2)
    end
    xlabel!("Time")
    ylabel!("Proportion")
    title!("Species Proportions Over Time")
end



end  # module SpeciesUtils
