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

using LinearAlgebra:    ×
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
    bField ::Array{wpFloat,4} # The magnetic field
    eField ::Array{wpFloat,4} # The eletric field
    ∇B     ::Array{wpFloat,4} # The gradient of the magnetic field
    ∇b̂     ::Array{wpFloat,5} # The gradient of the magnetic field
    ∇ExB   ::Array{wpFloat,5} # The gradient of the magnetic field
    xCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    yCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    zCoords::Vector{wpFloat} # The cartesian x-coordinates of the grid points
    domain ::Matrix{wpFloat} # Contains the extent of the numerical domain
    numdims::wpInt

    # Constructors
    #--------------------------------------------------------------------------
    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4},
                  ∇B     ::Array{wpFloat, 4},
                  ∇b̂     ::Array{wpFloat, 5},
                  ∇ExB   ::Array{wpFloat, 5},
                  xCoords::Vector{wpFloat},
                  yCoords::Vector{wpFloat},
                  zCoords::Vector{wpFloat}
                  )
        domain = [xCoords[1] xCoords[end]
                  yCoords[1] yCoords[end]
                  zCoords[1] zCoords[end]]
        numdims = 3
        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
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
        ∇B, ∇b̂, ∇ExB = compute∇s(bField, eField,
                                 xCoords,
                                 yCoords,
                                 zCoords)
        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
                   xCoords, yCoords, zCoords, 
                   domain, numdims)
    end # constructor

    function Mesh(bField ::Array{wpFloat, 4},
                  eField ::Array{wpFloat, 4}
                  )
        xCoords = collect(LinRange(0,1, size(bField)[2]))
        yCoords = collect(LinRange(0,1, size(bField)[3]))
        zCoords = collect(LinRange(0,1, size(bField)[4]))
        domain = [xCoords[1] xCoords[end]
                  yCoords[1] yCoords[end]
                  zCoords[1] zCoords[end]]
#        ∇B = compute∇B(bField, 
#                       xCoords,
#                       yCoords,
#                       zCoords)
        ∇B, ∇b̂, ∇ExB = compute∇s(bField, eField,
                                 xCoords,
                                 yCoords,
                                 zCoords)
        numdims = 3
        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
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
        #∇B = ∇(bField, dx, dy, derivate4thOrder) # Won't work. need norm3
        ∇B, ∇b̂, ∇ExB = compute∇s(bField, eField, # Won't work. Need norm3 etc.
                                 xCoords,
                                 yCoords,
                                 zCoords)
        numdims = 2
        # Won't work. missing ∇b̂ and ∇ExB
        return new(bField, eField, ∇B, xCoords, yCoords, numdims)
    end # constructor 

    function Mesh(bField ::Array{wpFloat, 2},
                  eField ::Array{wpFloat, 2}
                  )
        xCoords = LinRange(0,1, size(bField)[2])
        dx = xCoords[2] - xCoords[1]
        ∇B = derivateCentral(bField, dx, (1,0,0)) # Won't work. need norm2
        numdims = 1
        # Won't work. missing ∇b̂ and ∇ExB
        return new(bField, eField, ∇B, xCoords, numdims)
    end # constructor 
end # Struc

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


function compute∇s(bField ::Array{wpFloat, 4},
                   eField ::Array{wpFloat, 4},
                   xCoords::Vector{wpFloat},
                   yCoords::Vector{wpFloat},
                   zCoords::Vector{wpFloat}
                   )
    _, nx, ny, nz = size(bField)
    dx = xCoords[2] - xCoords[1]
    dy = yCoords[2] - yCoords[1]
    dz = zCoords[2] - zCoords[1]
    BB = norm4(bField)
    b̂ = zeros(size(bField))
    ExBdrift = zeros(size(bField))
    for i = 1:nx
        for j= 1:ny
            for k = 1:nz
                B⃗ = bField[:,i,j,k]
                E⃗ = eField[:,i,j,k]
                B = BB[i,j,k]
                b̂[:,i,j,k]  .= B⃗ ./ B
                ExBdrift[:,i,j,k] .= (E⃗ × B⃗) ./ B^2
            end
        end
    end
    ∇B = ∇(BB,  dx, dy, dz, derivate4thOrder) 
    ∇b̂ = ∇(b̂,  dx, dy, dz, derivate4thOrder)
    ∇ExBdrift= ∇(ExBdrift, dx, dy, dz, derivate4thOrder)
    return ∇B, ∇b̂, ∇ExBdrift
end # function compute∇s
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
