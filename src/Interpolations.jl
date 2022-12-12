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
    # Interpolate
    B = zeros(numDims) # The interpolated magnetic field
    E = zeros(numDims) # The interpolated electric field
    for l = 1:numDims # loop over dimensions
        B[l] = interpolator(mesh.bField[l, :, :, :], (i,j,k), (x,y,z))
        E[l] = interpolator(mesh.eField[l, :, :, :], (i,j,k), (x,y,z))
    end # loop over dimensions
    return B, E
    
end # function grid

end # module Interpolations
