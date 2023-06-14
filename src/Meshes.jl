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
using Random:           MersenneTwister

using Bifrost

using Schemes:          derivate4thOrder, derivateCentral, ∇, derivateUpwind
using Utilities:        norm4
using Constants
using Schemes


# Simple mesh
export Mesh
export amplifyBfield! # Amplifies the magnetic field of the mesh by a factor
export amplifyEfield! # Amplifies the elctric field of the mesh by a factor

#-------------#
# Structs     #
#-------------#-----------------------------------------------------------------
struct Mesh
    bField ::Array{T, 4} where {T<:Real} # The magnetic field
    eField ::Array{T, 4} where {T<:Real} # The eletric field
    ∇B     ::Array{T, 4} where {T<:Real} # The gradient of the magnetic field
    ∇b̂     ::Array{T, 5} where {T<:Real} # The gradient of the magnetic field
    ∇ExB   ::Array{T, 5} where {T<:Real} # The gradient of the magnetic field
    xCoords::Vector{T} where {T<:Real} # The cartesian x-coordinates of the grid points
    yCoords::Vector{T} where {T<:Real} # The cartesian x-coordinates of the grid points
    zCoords::Vector{T} where {T<:Real} # The cartesian x-coordinates of the grid points
    domain ::Matrix{T} where {T<:Real} # Contains the extent of the numerical domain
    numdims::Integer

    # Constructors
    #--------------------------------------------------------------------------
    function Mesh(bField ::Array{T, 4} where {T<:Real},
                  eField ::Array{T, 4} where {T<:Real},
                  ∇B     ::Array{T, 4} where {T<:Real},
                  ∇b̂     ::Array{T, 5} where {T<:Real},
                  ∇ExB   ::Array{T, 5} where {T<:Real},
                  xCoords::Vector{T} where {T<:Real},
                  yCoords::Vector{T} where {T<:Real},
                  zCoords::Vector{T} where {T<:Real}
                  )
        domain = [xCoords[1] xCoords[end]
                  yCoords[1] yCoords[end]
                  zCoords[1] zCoords[end]]
        numdims = 3
        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
                   xCoords, yCoords, zCoords, 
                   domain, numdims)
    end # constructor

    function Mesh(bField ::Array{T, 4} where {T<:Real},
                  eField ::Array{T, 4} where {T<:Real},
                  xCoords::Vector{T} where {T<:Real},
                  yCoords::Vector{T} where {T<:Real},
                  zCoords::Vector{T} where {T<:Real}
                  )
        domain = [xCoords[1] xCoords[length(xCoords)]
                  yCoords[1] yCoords[length(yCoords)]
                  zCoords[1] zCoords[length(zCoords)]]
        numdims = 3
        ∇B, ∇b̂, ∇ExB = compute∇s(bField, eField,
                                 xCoords,
                                 yCoords,
                                 zCoords,
                                 derivateUpwind
                                 )
        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
                   xCoords, yCoords, zCoords, 
                   domain, numdims)
    end # constructor

    function Mesh(bField ::Array{T, 4} where {T<:Real},
                  eField ::Array{T, 4} where {T<:Real}
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

    function Mesh(bField ::Array{T, 3} where {T<:Real},
                  eField ::Array{T, 3} where {T<:Real}
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

    function Mesh(bField ::Array{T, 2} where {T<:Real},
                  eField ::Array{T, 2} where {T<:Real}
                  )
        xCoords = LinRange(0,1, size(bField)[2])
        dx = xCoords[2] - xCoords[1]
        ∇B = derivateCentral(bField, dx, (1,0,0)) # Won't work. need norm2
        numdims = 1
        # Won't work. missing ∇b̂ and ∇ExB
        return new(bField, eField, ∇B, xCoords, numdims)
    end # constructor 

    """
         Using Bifrost input
    """
    function Mesh(
        expname::String,
        expdir ::String,
        snap   ::Integer,
        )

        numdims = 3
        #-----------------------------------------------------------------------
        # Load Bifrost-snap

        # Parse filenames
        snapstr = "_$snap"
        basename = string(expdir, "/", expname, "_$snap")
        idl_filename = string(basename, ".idl")
        snap_filename = string(basename, ".snap")
        aux_filename = string(basename, ".aux")
        mesh_filename = string(expdir, "/", expname, ".mesh")

        br_mesh = BifrostMesh(mesh_filename)
        br_params = br_read_params(idl_filename)
        br_snap = br_load_snapdata(snap_filename, br_params)
        br_aux = br_load_auxdata(aux_filename, br_params)

        auxvars = split(br_params["aux"])
        aux_avail = [0,                # η_total
                     0, 0, 0,  # ηx, ηy, ηz
                     0, 0, 0,  # Jx, Jy, Jz
                     0, 0, 0,  # Ex, Ey, Ez
                     ]
        for i = 1:length(auxvars)
            if auxvars[i] == "eta_total"
                aux_avail[1] = i
            elseif auxvars[i] == "etax"
                aux_avail[3] = i
            elseif auxvars[i] == "etay"
                aux_avail[4] = i
            elseif auxvars[i] == "etaz"
                aux_avail[2] = i
            elseif auxvars[i] == "ix"
                aux_avail[5] = i
            elseif auxvars[i] == "iy"
                aux_avail[6] = i
            elseif auxvars[i] == "iz"
                aux_avail[7] = i
            elseif auxvars[i] == "ex"
                aux_avail[8] = i
            elseif auxvars[i] == "ey"
                aux_avail[9] = i
            elseif auxvars[i] == "ez"
                aux_avail[10] = i
            end
        end

        #-----------------------------------------------------------------------
        # Allocate memory for simple variables
        meshsize = (br_mesh.mx, br_mesh.my, br_mesh.mz)
        wp_snap = Float32

        # Density
        ρ  = br_snap[:,:,:,1]
        # Bulk momentum
        px = br_snap[:,:,:,2]
        py = br_snap[:,:,:,3]
        pz = br_snap[:,:,:,4]
        # Magnetic field
        bx = br_snap[:,:,:,6]
        by = br_snap[:,:,:,7]
        bz = br_snap[:,:,:,8]

        dx = diff(br_mesh.x)
        dy = diff(br_mesh.y)
        dz = diff(br_mesh.z)

        domain = [br_mesh.x[1] br_mesh.x[end]
                  br_mesh.y[1] br_mesh.y[end]
                  br_mesh.z[1] br_mesh.z[end]]

        pbc_x = Bool(br_params["periodic_x"])
        pbc_y = Bool(br_params["periodic_y"])
        pbc_z = Bool(br_params["periodic_z"])

        #-----------------------------------------------------------------------
        # Destagger and scale variables to SI-units
        #
        # Vector quiantities are defined at the faces, while scalar quantities
        # are defined at cell centres. To get the bulk velcity at the centre
        # one thus have to interpolate the momentum to the cell centres.
        #
        # br_xup moves values half a grid in the x-direction.
        # params["u_u"] scales velocity from model/code-units to CGS units. 
        # params["u_b"] scales magnetic field from model/code-units to CGS units. 
        # cgs2SI_u scales velocity from CGS-units to SI-units
        # cgs2SI_b scales magnetic field from CGS-units to SI-units
        code2cgs_u = wp_snap(br_params["u_u"])
        code2cgs_B = wp_snap(br_params["u_B"])

        ux_cgs = code2cgs_u * br_xup(px, pbc_x) ./ ρ
        uy_cgs = code2cgs_u * br_yup(py, pbc_y) ./ ρ
        uz_cgs = code2cgs_u * br_zup(pz, pbc_z) ./ ρ
        ux_SI = ux_cgs * Constants.cgs2SI_u
        uy_SI = uy_cgs * Constants.cgs2SI_u
        uz_SI = uz_cgs * Constants.cgs2SI_u

        # De-stagger and scale magnetic field
        # bx_SI = bx_cgs * Constants.cgs2SI_b
        # by_SI = by_cgs * Constants.cgs2SI_b
        # bz_SI = bz_cgs * Constants.cgs2SI_b
        # Only scale magnetic field
        bx_SI = code2cgs_B * Constants.cgs2SI_b * bx
        by_SI = code2cgs_B * Constants.cgs2SI_b * by
        bz_SI = code2cgs_B * Constants.cgs2SI_b * bz
        
        #-----------------------------------------------------------------------
        # Compute the electric field and current density
        eField  = zeros(wp_snap, 3, meshsize...)
        calcEfield = true
        resistiveMHD = true

        if (aux_avail[8] > 0) & (aux_avail[9] > 0) & (aux_avail[10] > 0)
            eField[1,:,:,:] = br_aux[:,:,:,aux_avail[8]]
            eField[2,:,:,:] = br_aux[:,:,:,aux_avail[9]]
            eField[3,:,:,:] = br_aux[:,:,:,aux_avail[10]]
            calcEfield = false
        elseif  (aux_avail[2] > 0) & (aux_avail[3] > 0) & (aux_avail[4] > 0)
            # Calculate η_total? or J using components?
            println(string("Error: Computing E-field using η-components not ",
                           "implemented"))
            return
        elseif aux_avail[1] > 0
            η = br_aux[:,:,:,aux_avail[1]]
            # Do nothing
        else
            println("Warning: No resisstance nor electric field in aux-vars.")
            println("         Proceeds assuming ideal MHD.")
            resistiveMHD = false
        end

        # Bring magnetic field to cell centres
        bField  = zeros(wp_snap, 3, meshsize...)
        bField[1,:,:,:] = br_xup(bx_SI, pbc_x)
        bField[2,:,:,:] = br_xup(by_SI, pbc_y)
        bField[3,:,:,:] = br_xup(bz_SI, pbc_z)

#        # Temporary treatment of 2D meshes. Necessary for trilinear interpolations.
#        # Idealy, you want a bilinear_xz interpolation function to avoid
#        # interpolation in the third dimension.
#        if length(br_mesh.y) == 1
#            br_mesh.y = [br_mesh.y[1], br_mesh.y[1]]
#        end

        if calcEfield
            # This allocation is not needed if cross() allows components and
            # not only the vectors.
            bulkvel = zeros(wp_snap, 3, meshsize...)
            bulkvel[1,:,:,:] = ux_SI
            bulkvel[2,:,:,:] = uy_SI
            bulkvel[3,:,:,:] = uz_SI
            if resistiveMHD
                # Calculate current density
                if (aux_avail[5] > 0) & (aux_avail[6] > 0) & (aux_avail[7] > 0)
                    # Current density available in aux-variables
                    J = zeros(wp_snap, 3, meshsize...)
                    J[1,:,:,:] = br_aux[:,:,:,aux_avail[5]]
                    J[2,:,:,:] = br_aux[:,:,:,aux_avail[6]]
                    J[3,:,:,:] = br_aux[:,:,:,aux_avail[7]]
                elseif aux_avail[1] > 0
                    # Calculate J using the curl of the magnetic field.
                    # Use magnetic field on cell faces and upwind scheme
                    J = Schemes.curl(bx_SI, by_SI, bz_SI,
                                     br_mesh.x, br_mesh.y, br_mesh.z,
                                     Schemes.derivateUpwind
                                     )/Constants.μ_0
                end
                eField = J - cross(bulkvel, bField)
                eField[1,:,:,:] .*= η
                eField[2,:,:,:] .*= η
                eField[3,:,:,:] .*= η
            else
                eField = -cross(bulkvel, bField)
            end
        end

        #-----------------------------------------------------------------------
        # Define/compute other parameters/variables.

        ∇B, ∇b̂, ∇ExB = compute∇s(bField, eField,
                                 br_mesh.x,
                                 br_mesh.y,
                                 br_mesh.z,
                                 derivateUpwind,
                                 )

        #-------------------------------------------------------

        return new(bField, eField, ∇B, ∇b̂, ∇ExB,
                   br_mesh.x, br_mesh.y, br_mesh.z,
                   domain, numdims)
    end # constructor Mesh

end # Struc

#-------------------#
# Utility functions #
#-------------------------------------------------------------------------------
function compute∇B(bField ::Array{T, 4} where {T<:Real},
                   xCoords::Vector{T} where {T<:Real},
                   yCoords::Vector{T} where {T<:Real},
                   zCoords::Vector{T} where {T<:Real}
                   )
    dx = xCoords[2] - xCoords[1]
    dy = yCoords[2] - yCoords[1]
    dz = zCoords[2] - zCoords[1]
    ∇B = ∇(norm4(bField), dx, dy, dz, derivate4thOrder) 
    return ∇B
end # function compute∇B


function compute∇s(bField ::Array{T, 4} where {T<:Real},
                   eField ::Array{T, 4} where {T<:Real},
                   xCoords::Vector{T} where {T<:Real},
                   yCoords::Vector{T} where {T<:Real},
                   zCoords::Vector{T} where {T<:Real}
                   )
    wfp = typeof(bField[1])
    _, nx, ny, nz = size(bField)
    dx = xCoords[2] - xCoords[1]
    dy = yCoords[2] - yCoords[1]
    dz = zCoords[2] - zCoords[1]
    BB = norm4(bField)
    b̂ = zeros(wfp, 3, nx, ny, nz)
    ExBdrift = zeros(wfp, 3, nx, ny, nz)
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
    ∇B = ∇(BB, dx, dy, dz, derivate4thOrder) 
    ∇b̂ = ∇(b̂,  dx, dy, dz, derivate4thOrder)
    ∇ExBdrift= ∇(ExBdrift, dx, dy, dz, derivate4thOrder)
    return ∇B, ∇b̂, ∇ExBdrift
end # function compute∇s
#|
function compute∇s(bField ::Array{T, 4} where {T<:Real},
                   eField ::Array{T, 4} where {T<:Real},
                   xCoords::Vector{T} where {T<:Real},
                   yCoords::Vector{T} where {T<:Real},
                   zCoords::Vector{T} where {T<:Real},
                   scheme ::Function
                   )
    wfp = typeof(bField[1])
    _, nx, ny, nz = size(bField)
    BB = norm4(bField)
    b̂ = zeros(wfp, 3, nx, ny, nz)
    ExBdrift = zeros(wfp, 3, nx, ny, nz)
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
    ∇B = ∇(BB, xCoords, yCoords, zCoords, scheme) 
    ∇b̂ = ∇(b̂,  xCoords, yCoords, zCoords, scheme)
    ∇ExBdrift= ∇(ExBdrift, xCoords, yCoords, zCoords, scheme)
    return ∇B, ∇b̂, ∇ExBdrift
end # function compute∇s
#------------------#
# Mesh set-methods #
#-------------------------------------------------------------------------------
function amplifyBfield!(mesh  ::Mesh,
                        factor::Real)
    @. mesh.bField *= factor
    mesh.∇B = compute∇B(mesh.bField, 
                        mesh.xCoords,
                        mesh.yCoords,
                        mesh.zCoords)
end # function amplifyBfield
                   

function amplifyEfield!(mesh  ::Mesh,
                        factor::Real)
    @. mesh.eField *= factor
end # function amplifyEfield


end # module Meshes
