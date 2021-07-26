export interaction_matrix, PowerLaw

abstract type Interaction end

"""
    interaction_matrix(interaction, geometry, positions)
    interaction_matrix(interaction, distances)

Compute the interaction matrix from given distances or geometry and positions.
The latter is only supported for anistrope interactions as these do only depend on the distance.
"""
function interaction_matrix end


abstract type AnisotropeInteraction <: Interaction end

interaction_matrix(interaction::AnisotropeInteraction, geometry::Geometry, positions) = interaction_matrix(interaction, distance_matrix(geometry, positions))

interaction_matrix(interaction::AnisotropeInteraction, distances) = _interaction_strength.(Ref(interaction), distances)

"""
    _interaction_strength(<:AnistropeInteraction, distance)
Compute interaction strength given the distance.
"""
function _interaction_strength end

"""
    PowerLaw{T}

Interaction strength scales proportional to distance to the power -α.

# Fields
- `α::T`: Exponent
"""
struct PowerLaw <: AnisotropeInteraction
    α::Float64
    PowerLaw(α) = new(convert(Float64, α))
end

_interaction_strength(int::PowerLaw, distance::Float64) = float(distance > 0 ? distance^-int.α : 0)