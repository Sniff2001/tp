#-------------------------------------------------------------------------------
# Created 02.12.22
# email: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                Meshes.jl
#
#-------------------------------------------------------------------------------
# Module containing mesh-structs and methods
#-------------------------------------------------------------------------------

module Meshes

using WorkingPrecision: wpFloat, wpInt
using Schemes:          derivate4thOrder, derivateCentral, ∇, norm4

# Simple mesh
export Mesh

#-------------#
# Structs     #
#-------------#-----------------------------------------------------------------
struct Mesh
    bField ::AbstractArray{wpFloat}   # The magnetic field
    eField ::AbstractArray{wpFloat}   # The eletric field
    ∇B     ::Vector{Array{wpFloat, 3}} # The gradient of the magnetic field
    xCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    yCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    zCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points

    # Constructors
    #--------------------------------------------------------------------------
    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4},
                  ∇B     ::Array{wpFloat, 3},
                  xCoords::Vector{wpFloat},
                  yCoords::Vector{wpFloat},
                  zCoords::Vector{wpFloat}
                  )
        return new(bField, eField, ∇B, xCoords, yCoords, zCoords)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4},
                  xCoords::Vector{wpFloat},
                  yCoords::Vector{wpFloat},
                  zCoords::Vector{wpFloat}
                  )
        dx = xCoords[2] - xCoords[1]
        dy = yCoords[2] - yCoords[1]
        dz = zCoords[2] - zCoords[1]
        ∇B = ∇(norm4(bField), dx, dy, dz, derivate4thOrder) 
        return new(bField, eField, ∇B, xCoords, yCoords, zCoords)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4}
                  )
        xCoords = LinRange(0,1, size(bField)[2])
        yCoords = LinRange(0,1, size(bField)[3])
        zCoords = LinRange(0,1, size(bField)[4])
        dx = xCoords[2] - xCoords[1]
        dy = yCoords[2] - yCoords[1]
        dz = zCoords[2] - zCoords[1]
        ∇B = ∇(norm4(bField), dx, dy, dz, derivate4thOrder) 
        return new(bField, eField, ∇B, xCoords, yCoords, zCoords)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 3},
                  eField ::Array{wpFloat, 3}
                  )
        xCoords = LinRange(0,1, size(bField)[2])
        yCoords = LinRange(0,1, size(bField)[3])
        dx = xCoords[2] - xCoords[1]
        dy = yCoords[2] - yCoords[1]
        ∇B = ∇(bField, dx, dy, derivate4thOrder) # Won't work. need norm3
        return new(bField, eField, ∇B, xCoords, yCoords)
    end # constructor 

    function Mesh(bField ::Array{wpFloat, 2},
                  eField ::Array{wpFloat, 2}
                  )
        xCoords = LinRange(0,1, size(bField)[2])
        dx = xCoords[2] - xCoords[1]
        ∇B = derivateCentral(bField, dx, (1,0,0)) # Won't work. need norm2
        return new(bField, eField, ∇B, xCoords)
    end # constructor 
end # Struct
#-------------------------------------------------------------------------------

end # Module
