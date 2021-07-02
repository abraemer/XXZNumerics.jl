export interaction_matrix, PowerLaw

abstract type Interaction end

"""
    interaction_matrix(interaction, geometry, positions)
    interaction_matrix(interaction, distances)

Compute the interaction matrix from given distances or geometry and positions.
"""
function interaction_matrix(::Interaction, ::Geometry, positions) end


abstract type AnisotropeInteraction <: Interaction end

interaction_matrix(interaction::AnisotropeInteraction, geometry::Geometry, positions) = interaction_matrix(interaction, distance_matrix(geometry, positions))

interaction_matrix(interaction::AnisotropeInteraction, distances) = _interaction_strength.(Ref(interaction), distances)

"""
    PowerLaw{T}

Interaction strength scales proportional to distance to the power -α.

# Fields
- `α::T`: Exponent
"""
struct PowerLaw{T} <: AnisotropeInteraction
    α::T
end

_interaction_strength(int::PowerLaw, distance) = distance > 0 ? distance^-int.α : 0