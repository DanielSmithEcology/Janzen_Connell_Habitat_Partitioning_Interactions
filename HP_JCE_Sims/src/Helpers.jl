# src/Helpers.jl
module Helpers

export closest_greater_index

"""
    closest_greater_index(values::Vector{T}, target::T)

Return the index of the value in `values` that is the smallest value strictly greater than `target`.
Returns 0 if no such value exists.
"""
function closest_greater_index(values::Vector{T}, target::T) where T
    closest_idx = 0
    min_diff = typemax(T)
    for (i, val) in pairs(values)
        if val > target && (val - target) < min_diff
            closest_idx = i
            min_diff = val - target
        end
    end
    return closest_idx
end

end  # module
