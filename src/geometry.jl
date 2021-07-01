export sample_blockaded, distance, distance_matrix, Box, BoxPBC, NoisyChain, NoisyChainPBC

## TODO: abstract Lattice geometries and put chains as subtypes
##       lattices need to know at which site to generate a new point
## TODO: orthogonalize PBC types

abstract type Geometry end

"""
    generate_point(geometry)

Generate a new (random) point in the geometry.
"""
function generate_point(geometry::Geometry)
    error("Not implemented!")
end

_euclidean_point(lengths) = rand(length(lengths)) .* lengths

"""
    distance(geometry, p1, p2)

Compute the distance of two points p1 and p2 within the given geometry.
"""
function distance(geometry::Geometry, p1, p2)
    error("Not implemented!")
end

## two default implementations

function _euclidean(p1, p2)
    sqrt(sum(@. (p1 - p2)^2 ))
end

function _euclidean_pbc(p1, p2, lengths)
    diff = @. abs(mod(p1, lengths) - mod(p2, lengths))
    dim_mask = diff .> lengths./2
    diff[dim_mask] .-= lengths[dim_mask]
    sqrt(sum(diff.^2))
end

"""
    distance_matrix(geometry, points)

Compute the distance pairwise between all points.
"""
function distance_matrix(geometry::Geometry, points)
    N = length(points)
    res = zeros(Float64, (N, N))
    for i in 1:N
        for j in i+1:N
            res[i,j] = distance(geometry, points[i], points[j])
        end
    end
    Symmetric(res)
end

"""
    sample_blockaded(geometry, N; blockade_radius, max_iter)

Generate N points with distance blockade_radius (default 1) within the given geometry.
"""
function sample_blockaded(geometry, N; blockade_radius=1, max_iter=1000)
    res = [generate_point(geometry)]
    sizehint!(res, N)
    while length(res) < N
        push!(res, _find_new_point(geometry, res, blockade_radius, max_iter))
    end
    res
end

function _find_new_point(geometry, points, blockade_radius, max_iter)
    for _ in 1:max_iter
        p = generate_point(geometry)
        if all(map(p2 -> distance(geometry, p, p2) > blockade_radius, points))
            return p
        end
    end
    error("Iteration limit reached! Could not find a $(length(points)+1)th point in $geometry geometry with $max_iter iterations.")
end

"""
    Box{T}

Geometry representing a simple box with given lengths.

# Fields
- `lengths::Vector{T}`: stores the lengths of the individual dimensions
"""
struct Box{T} <: Geometry 
    lengths::Vector{T}
end

generate_point(box::Box) = _euclidean_point(box.lengths)
distance(::Box, p1, p2) = _euclidean(p1, p2)


"""
    BoxPBC{T}

Geometry representing a box with given lengths respecting periodic boundary conditions.

# Fields
- `lengths::Vector{T}`: stores the lengths of the individual dimensions
"""
struct BoxPBC{T} <: Geometry
    lengths::Vector{T}
end

distance(box::BoxPBC, p1, p2) = _euclidean_pbc(p1, p2, box.lengths)
generate_point(box::BoxPBC) = _euclidean_point(box.lengths)

"""
    NoisyChain

Represents a 1D lattice with given spacing and length.
Each position is subjected to noise drawn from a box distribution on [-σ, σ].

# Fields
- `L::Int64`: number of sites in the chain
- `spacing::Float64`: (mean) distance between two sites
- `σ::Float64`: width of box distribution for position noise
"""
struct NoisyChain
    L::Int64
    spacing::Float64
    σ::Float64
end

distance(::NoisyChain, p1, p2) = _euclidean(p1, p2)

generate_point(chain::NoisyChain, site=1) = chain.σ * (2*rand()-1) + site*chain.spacing

"""
    NoisyChainPBC

Represents a 1D periodic lattice with given spacing and length.
Each position is subjected to noise drawn from a box distribution on [-σ, σ].

# Fields
- `L::Int64`: number of sites in the chain
- `spacing::Float64`: (mean) distance between two sites
- `σ::Float64`: width of box distribution for position noise
"""
struct NoisyChainPBC
    L::Int64
    spacing::Float64
    σ::Float64
end

_chain_length(chain::NoisyChainPBC) = chain.L * chain.spacing

distance(chain::NoisyChainPBC, p1, p2) = _euclidean_pbc(p1, p2, _chain_length(chain))

generate_point(chain::NoisyChainPBC, site=1) = mod(chain.σ * (2*rand()-1) + site*chain.spacing,  _chain_length(chain))

## chain types need to now at which site to generate the points
function _find_new_point(chain::Union{NoisyChain, NoisyChainPBC}, points, blockade_radius, max_iter)
    site = length(points) + 1
    for _ in 1:max_iter
        p = generate_point(chain, site)
        if all(map(p2 -> distance(chain, p, p2) > blockade_radius, points))
            return p
        end
    end
    error("Iteration limit reached! Could not find a $(length(points)+1)th point in $chain geometry with $max_iter iterations.")
end
