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

using PyPlot

using WorkingPrecision: wpFloat, wpInt
using Schemes:          derivate4thOrder, derivateCentral, ∇, norm4

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
        numdims = 4
        return new(bField, eField, ∇B, xCoords, yCoords, zCoords, numdims)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4},
                  xCoords::Vector{wpFloat},
                  yCoords::Vector{wpFloat},
                  zCoords::Vector{wpFloat}
                  )
        ∇B = compute∇B(bField, 
                       xCoords,
                       yCoords,
                       zCoords)
        return new(bField, eField, ∇B, xCoords, yCoords, zCoords, 4)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4}
                  )
        xCoords = collect(LinRange(0,1, size(bField)[2]))
        yCoords = collect(LinRange(0,1, size(bField)[3]))
        zCoords = collect(LinRange(0,1, size(bField)[4]))
        ∇B = compute∇B(bField, 
                       xCoords,
                       yCoords,
                       zCoords)
        numdims = 4
        return new(bField, eField, ∇B, xCoords, yCoords, zCoords, numdims)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 3},
                  eField ::Array{wpFloat, 3}
                  )
        xCoords = collect(LinRange(0,1, size(bField)[2]))
        yCoords = collect(LinRange(0,1, size(bField)[3]))
        dx = xCoords[2] - xCoords[1]
        dy = yCoords[2] - yCoords[1]
        ∇B = ∇(bField, dx, dy, derivate4thOrder) # Won't work. need norm3
        numdims = 3
        return new(bField, eField, ∇B, xCoords, yCoords, numdims)
    end # constructor 

    function Mesh(bField ::Array{wpFloat, 2},
                  eField ::Array{wpFloat, 2}
                  )
        xCoords = LinRange(0,1, size(bField)[2])
        dx = xCoords[2] - xCoords[1]
        ∇B = derivateCentral(bField, dx, (1,0,0)) # Won't work. need norm2
        numdims = 2
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


#---------------#
# Mesh plotting #
#-------------------------------------------------------------------------------
function quiverslice(mesh::Mesh, 
                     normal::String, 
                     point::wpInt,
                     field::String="B")
    if field == "B"
        f = mesh.bField
    elseif field == "E"
        f = mesh.eField
    else
        println("Error: Invalid field.")
    end
    if normal == "x"
        quiver(mesh.yCoords,
               mesh.zCoords,
               transpose(f[2, point, :, :]),
               transpose(f[3, point, :, :])
               )
        xlabel("y")
        ylabel("z")
    elseif normal == "y"
        quiver(mesh.xCoords,
               mesh.zCoords,
               transpose(f[1, :, point, :]),
               transpose(f[3, :, point, :])  
               )
        xlabel("x")
        ylabel("z")
    elseif normal == "z"
        quiver(mesh.xCoords,
               mesh.yCoords,
               transpose(f[1, :, :, point]),
               transpose(f[2, :, :, point])  
               )
        xlabel("x")
        ylabel("y")
    else
        println("Error: Plane not valid.")
    end
end # function quiverslice


function streamplotslice(mesh::Mesh, 
                         normal::String, 
                         point::wpInt,
                         field::String="B")
    if field == "B"
        f = mesh.bField
    elseif field == "E"
        f = mesh.eField
    else
        println("Error: Invalid field.")
    end
    if normal == "x"
        streamplot(mesh.yCoords,
               mesh.zCoords,
               transpose(f[2, point, :, :]),
               transpose(f[3, point, :, :])
               )
        xlabel("y")
        ylabel("z")
    elseif normal == "y"
        streamplot(mesh.xCoords,
               mesh.zCoords,
               transpose(f[1, :, point, :]),
               transpose(f[3, :, point, :])  
               )
        xlabel("x")
        ylabel("z")
    elseif normal == "z"
        streamplot(mesh.xCoords,
               mesh.yCoords,
               transpose(f[1, :, :, point]),
               transpose(f[2, :, :, point])  
               )
        xlabel("x")
        ylabel("y")
    else
        println("Error: Plane not valid.")
    end
end # function streamplotslice

end # module Meshes
