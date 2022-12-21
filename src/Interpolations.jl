#-------------------------------------------------------------------------------
#
# Created 12.12.22
# Author: e.s.oyre@astro.uio.no
#
# Last edited: 13.12.22
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

export locateCell
export trilinear

function grid(mesh       ::Mesh,
              interpolator::Function,
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
        mid  = floor(wpInt, (high - low)/2) + low
        if point > coords[mid]
            low = mid
        else
            high = mid
        end # if
    end # while
    return low 
end # function locateCell


function trilinear(mesh   ::Mesh, 
                   (i,j,k)::Tuple{wpInt, wpInt, wpInt},
                   (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
                   )

    t = (x - mesh.xCoords[i])/(mesh.xCoords[i + 1] - mesh.xCoords[i])
    u = (y - mesh.yCoords[j])/(mesh.yCoords[j + 1] - mesh.yCoords[j])
    v = (z - mesh.zCoords[k])/(mesh.zCoords[k + 1] - mesh.zCoords[k])

    B0 = mesh.bField[:,   i,   j,   k]
    B1 = mesh.bField[:, i+1,   j,   k]
    B2 = mesh.bField[:, i+1, j+1,   k]
    B3 = mesh.bField[:,   i, j+1,   k]
    B4 = mesh.bField[:,   i,   j, k+1]
    B5 = mesh.bField[:, i+1,   j, k+1]
    B6 = mesh.bField[:, i+1, j+1, k+1]
    B7 = mesh.bField[:,   i, j+1, k+1]

    E0 = mesh.eField[:,   i,   j,   k]
    E1 = mesh.eField[:, i+1,   j,   k]
    E2 = mesh.eField[:, i+1, j+1,   k]
    E3 = mesh.eField[:,   i, j+1,   k]
    E4 = mesh.eField[:,   i,   j, k+1]
    E5 = mesh.eField[:, i+1,   j, k+1]
    E6 = mesh.eField[:, i+1, j+1, k+1]
    E7 = mesh.eField[:,   i, j+1, k+1]

    f0 = (1 - t)*(1 - u)*(1 - v)
    f1 = t*(1 - u)*(1 - v)
    f2 = t*u*(1 - v)
    f3 = (1 - t)*u*(1 - v)
    f4 = (1 - t)*(1 - u)*v
    f5 = t*(1 - u)*v
    f6 = t*u*v
    f7 = (1 - t)*u*v

    B = f0*B0 + f1*B1 + f2*B2 + f3*B3 + f4*B4 + f5*B5 + f6*B6 + f7*B7
    E = f0*E0 + f1*E1 + f2*E2 + f3*E3 + f4*E4 + f5*E5 + f6*E6 + f7*E7

    return B, E

end # function trinlinear

    
end # module Interpolations
