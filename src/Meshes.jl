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
using Schemes:          derivate4thOrder, derivateCentral, ∇
using Utilities:        norm4

# Simple mesh
export Mesh
export amplifyBfield! # Amplifies the magnetic field of the mesh by a factor
export amplifyEfield! # Amplifies the elctric field of the mesh by a factor

#-------------#
# Structs     #
#-------------#-----------------------------------------------------------------
struct Mesh
    bField ::AbstractArray{wpFloat} # The magnetic field
    eField ::AbstractArray{wpFloat} # The eletric field
    ∇B     ::AbstractArray{wpFloat} # The gradient of the magnetic field
    xCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    yCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    zCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    domain ::Matrix{wpFloat} # Contains the extent of the numerical domain
    numdims::wpInt

    # Constructors
    #--------------------------------------------------------------------------
    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4},
                  ∇B     ::Array{wpFloat, 3},
                  xCoords::Vector{wpFloat},
                  yCoords::Vector{wpFloat},
                  zCoords::Vector{wpFloat}
                  )
        domain = [xCoords[1] xCoords[length(xCoords)]
                  yCoords[1] yCoords[length(yCoords)]
                  zCoords[1] zCoords[length(zCoords)]]
        numdims = 3
        return new(bField, eField, ∇B, 
                   xCoords, yCoords, zCoords, 
                   domain, numdims)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4},
                  xCoords::Vector{wpFloat},
                  yCoords::Vector{wpFloat},
                  zCoords::Vector{wpFloat}
                  )
        domain = [xCoords[1] xCoords[length(xCoords)]
                  yCoords[1] yCoords[length(yCoords)]
                  zCoords[1] zCoords[length(zCoords)]]
        numdims = 3
        ∇B = compute∇B(bField, 
                       xCoords,
                       yCoords,
                       zCoords)
        return new(bField, eField, ∇B, 
                   xCoords, yCoords, zCoords, 
                   domain, numdims)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4}
                  )
        xCoords = collect(LinRange(0,1, size(bField)[2]))
        yCoords = collect(LinRange(0,1, size(bField)[3]))
        zCoords = collect(LinRange(0,1, size(bField)[4]))
        domain = [xCoords[1] xCoords[length(xCoords)]
                  yCoords[1] yCoords[length(yCoords)]
                  zCoords[1] zCoords[length(zCoords)]]
        ∇B = compute∇B(bField, 
                       xCoords,
                       yCoords,
                       zCoords)
        numdims = 3
        return new(bField, eField, ∇B, 
                   xCoords, yCoords, zCoords, 
                   domain, numdims)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 3},
                  eField ::Array{wpFloat, 3}
                  )
        xCoords = collect(LinRange(0,1, size(bField)[2]))
        yCoords = collect(LinRange(0,1, size(bField)[3]))
        dx = xCoords[2] - xCoords[1]
        dy = yCoords[2] - yCoords[1]
        ∇B = ∇(bField, dx, dy, derivate4thOrder) # Won't work. need norm3
        numdims = 2
        return new(bField, eField, ∇B, xCoords, yCoords, numdims)
    end # constructor 

    function Mesh(bField ::Array{wpFloat, 2},
                  eField ::Array{wpFloat, 2}
                  )
        xCoords = LinRange(0,1, size(bField)[2])
        dx = xCoords[2] - xCoords[1]
        ∇B = derivateCentral(bField, dx, (1,0,0)) # Won't work. need norm2
        numdims = 1
        return new(bField, eField, ∇B, xCoords, numdims)
    end # constructor 
end # Struct

#-------------------#
# Utility functions #
#-------------------------------------------------------------------------------
function compute∇B(bField ::Array{wpFloat, 4},
                   xCoords::Vector{wpFloat},
                   yCoords::Vector{wpFloat},
                   zCoords::Vector{wpFloat}
                   )
    dx = xCoords[2] - xCoords[1]
    dy = yCoords[2] - yCoords[1]
    dz = zCoords[2] - zCoords[1]
    ∇B = ∇(norm4(bField), dx, dy, dz, derivate4thOrder) 
    return ∇B
end # function compute∇B

#------------------#
# Mesh set-methods #
#-------------------------------------------------------------------------------
function amplifyBfield!(mesh  ::Mesh,
                        factor::wpFloat)
    @. mesh.bField *= factor
    mesh.∇B = compute∇B(mesh.bField, 
                        mesh.xCoords,
                        mesh.yCoords,
                        mesh.zCoords)
end # function amplifyBfield
                   

function amplifyEfield!(mesh  ::Mesh,
                        factor::wpFloat)
    @. mesh.eField *= factor
end # function amplifyEfield


end # module Meshes
