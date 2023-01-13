#-------------------------------------------------------------------------------
# Created 12.12.22
# Author: e.s.oyre@astro.uio.no
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

function grid(mesh        ::Mesh,
              interpolator::Function,
              pos         ::Vector{wpFloat}
              )
    numDims = 3
    x, y, z = pos
    # Find which cell the position is in, and store the indices of its corner.
    i = locateCell(mesh.xCoords, x)
    j = locateCell(mesh.yCoords, y)
    k = locateCell(mesh.zCoords, z)

    fields = interpolator(mesh, (i,j,k), (x,y,z))
    
    return fields, i, j, k

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
    while  high > low + 1 # While we haven't found the cell
        mid  = floor(wpInt, (high + low)/2)
        if point > coords[mid]
            low = mid 
        elseif point < coords[mid]
            high = mid
	else
	    low = mid - 1
	    high = mid
        end # if
    end # while
    return low 
end # function locateCell


function trilinear(mesh   ::Mesh,
                   (i,j,k)::Tuple{wpInt, wpInt, wpInt},
                   (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
                   )
    coefficients = trilinearCoefficients(mesh.xCoords,
                                         mesh.yCoords,
                                         mesh.xCoords,
                                         (i,j,k),
                                         (x,y,z))
    B = trilinearComputeField(mesh.bField,
                              (i,j,k),
                              coefficients)
    E = trilinearComputeField(mesh.eField,
                              (i,j,k),
                              coefficients)
    return B, E
end # function trilinear

function trilinearGCA(mesh   ::Mesh,
                      (i,j,k)::Tuple{wpInt, wpInt, wpInt},
                      (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
                      )
    coefficients = trilinearCoefficients(mesh.xCoords,
                                         mesh.yCoords,
                                         mesh.xCoords,
                                         (i,j,k),
                                         (x,y,z))
    B = trilinearComputeField(mesh.bField,
                              (i,j,k),
                              coefficients)
    E = trilinearComputeField(mesh.bField,
                              (i,j,k),
                              coefficients)
    ∇B = trilinearComputeField(mesh.∇B,
                              (i,j,k),
                              coefficients)
    return B, E, ∇B
end # function trilinearGCA

function trilinearCoefficients(xCoords::Vector{wpFloat},
                               yCoords::Vector{wpFloat},
                               zCoords::Vector{wpFloat},
                               (i,j,k)::Tuple{wpInt, wpInt, wpInt},
                               (x,y,z)::Tuple{wpFloat, wpFloat, wpFloat}
                               )
    t = (x - xCoords[i])/(xCoords[i + 1] - xCoords[i])
    u = (y - yCoords[j])/(yCoords[j + 1] - yCoords[j])
    v = (z - zCoords[k])/(zCoords[k + 1] - zCoords[k])
    #
    c0 = (1 - t)*(1 - u)*(1 - v)
    c1 = t*(1 - u)*(1 - v)
    c2 = t*u*(1 - v)
    c3 = (1 - t)*u*(1 - v)
    c4 = (1 - t)*(1 - u)*v
    c5 = t*(1 - u)*v
    c6 = t*u*v
    c7 = (1 - t)*u*v
    #
    return c0, c1, c2, c3, c4, c5, c6, c7
end # function trinlinearCoefficients

#_
function trilinearComputeField(field  ::Array{wpFloat, 4},
                               (i,j,k)::Tuple{wpInt, wpInt, wpInt},
                               c      ::NTuple{8, wpFloat}
                               )
    c0, c1, c2, c3, c4, c5, c6, c7 = c
    A0 = field[:,   i,   j,   k]
    A1 = field[:, i+1,   j,   k]
    A2 = field[:, i+1, j+1,   k]
    A3 = field[:,   i, j+1,   k]
    A4 = field[:,   i,   j, k+1]
    A5 = field[:, i+1,   j, k+1]
    A6 = field[:, i+1, j+1, k+1]
    A7 = field[:,   i, j+1, k+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 + c4*A4 + c5*A5 + c6*A6 + c7*A7
    return A
end # function trilinearComputeField

end # module Interpolations
