#-------------------------------------------------------------------------------
#
# Created 02.12.22
# By Eilif @RoCS UiO.
# email: e.s.oyre@astro.uio.no
#
# Last edited: 07.12.22
#
#-------------------------------------------------------------------------------
#
#                Meshes.jl
#
#-------------------------------------------------------------------------------
# Module containing mesh-structs and methods
#-------------------------------------------------------------------------------

module Meshes

using WorkingPrecision: wpFloat, wpInt

# Simple mesh
export Mesh

#-------------#
# Structs     #
#-------------#-----------------------------------------------------------------
struct Mesh
    bField ::AbstractArray{wpFloat} # The magnetic field
    eField ::AbstractArray{wpFloat} # The eletric field
    xCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    yCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    zCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points

    # Constructors
    #--------------------------------------------------------------------------
    function Mesh(bField ::Array{wpFloat, 2},
                  eField ::Array{wpFloat, 2})
        xCoords = LinRange(0,1, size(bField)[2])
        return new(bField, eField, xCoords)
    end # constructor 

    function Mesh(bField ::Array{wpFloat, 3},
                  eField ::Array{wpFloat, 3})
        xCoords = LinRange(0,1, size(bField)[2])
        yCoords = LinRange(0,1, size(bField)[3])
        return new(bField, eField, xCoords, yCoords)
    end # constructor 

    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4})
        xCoords = LinRange(0,1, size(bField)[2])
        yCoords = LinRange(0,1, size(bField)[3])
        zCoords = LinRange(0,1, size(bField)[4])
        return new(bField, eField, xCoords, yCoords, zCoords)
    end # constructor
end # Struct
#-------------------------------------------------------------------------------

end # Module
