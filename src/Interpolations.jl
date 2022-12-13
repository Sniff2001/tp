#-------------------------------------------------------------------------------
#
# Created 12.12.22
# Author: e.s.oyre@astro.uio.no
#
# Last edited: 12.12.22
#
#-------------------------------------------------------------------------------
#
#             Interpolations.jl
#
#-------------------------------------------------------------------------------
# Module containing interpolation methods.
#-------------------------------------------------------------------------------

module Interpolations

using WorkingPrecision: wpFloat, wpInt
using Meshes

function grid(mesh       ::Mesh,
              interpolator::Function
              pos        ::Vector{wpFloat}
              )
    numDims = 3
    x, y, z = pos
    # Find which cell the position is in, and store the indices of its corner.
    i = locateCell(mesh.xCoords, x)
    j = locateCell(mesh.yCoords, y)
    k = locateCell(mesh.zCoords, z)

    B, E = interpolator(mesh, (i,j,k), (x,y,z))
    
    return B, E

end # function grid


"""
    locateCell(coords, point)

Finds the position of `point` in the vector `coords` using binary search.
Returns its lower neighbour.
"""
function locateCell(coords::Vector{wpFloat},
                    point ::Number
                    )
    # Initial bounds
    low = 1
    high = length(coords)
    while  (high - low) > 1 # While we haven't found the cell
        mid  = convert(wpInt, (high - low)/2)
        if point > coords[mid]
            low = mid
        else
            high = mid
        end # if
    end # while
    return low 
end # function locateCell

    
end # module Interpolations
