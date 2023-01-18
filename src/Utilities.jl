#-------------------------------------------------------------------------------
# Created 17.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                Utilities.jl
#
#-------------------------------------------------------------------------------
# Module containing the utility functions.
#-------------------------------------------------------------------------------


module Utilities

# Standard libraries
using LinearAlgebra:    norm
# Internal libraries
using WorkingPrecision: wpFloat, wpInt


#----------------#
# Linear albebra #
#----------------#--------------------------------------------------------------
"""
    norm2(field, axis)
Calculates the p=2 norm of the vectors in a 1D vector field, i.e. the field
strength. The argument 'axis' determines whether the vector components are
stored in the first or second dimension of the array storing the field.
"""
function norm2(field::Matrix{wpFloat},
               axis ::wpInt=1
               )
    dims = size(field)
    if axis == 1
        fieldstrength = zeros(dims[2])
        for i = 1:dims[2]
            fieldstrength[i] = norm(field[:, i])
        end # loop i
    elseif axis == 2
        fieldstrength = zeros(dims[1])
        for i = 1:dims[1]
            fieldstrength[i] = norm(field[i, :])
        end
    else
        println("Error: Your axes are wierd...")
    end # if
    return fieldstrength
end # function norm2


"""
    norm3(field)
Calculates the p=2 norm of the vectors in a 2D vector field, i.e. the field
strength. The function assumes the vector components are store in the first
dimension of the field array.
"""
function norm3(field::Array{wpFloat, 3})
    dims = size(field)
    fieldstrength = zeros(dims[2:3])
    for i = 1:dims[2]
        for j = 1:dims[3]
            fieldstrength[i,j] = norm(field[:, i, j])
        end # loop j
    end # loop i
    return fieldstrength
end # function norm2


"""
    norm4(field, axis)
Calculates the p=2 norm of the vectors in a 3D vector field, i.e. the field
strength. The argument 'axis' determines whether the vector components are
stored in the first or fourth dimension of the array storing the field.
"""
function norm4(field::Array{wpFloat, 4},
               axis ::wpInt=1
               )
    if axis == 1
        dims = size(field[1,:,:,:])
        fieldstrength = zeros(dims)
        for i = 1:dims[1]
            for j = 1:dims[2]
                for k = 1:dims[3]
                    fieldstrength[i,j,k] = √(field[1,i,j,k]^2 + 
                                             field[2,i,j,k]^2 +
                                             field[3,i,j,k]^2)
                end # loop k
            end # look j
        end # loop i
    elseif axis == 4
        dims = size(field[:,:,:,1])
        fieldstrength = zeros(dims)
        for i = 1:dims[1]
            for j = 1:dims[2]
                for k = 1:dims[3]
                    fieldstrength[i,j,k] = √(field[i,j,k,1]^2 + 
                                             field[i,j,k,2]^2 +
                                             field[i,j,k,3]^2)
                end # loop k
            end # look j
        end # loop i
    else
        println("Error: Yours axes are wierd...")
    end # if
    return fieldstrength
end # function norm4


#----------------#
# Distributions  #
#----------------#--------------------------------------------------------------
function normaldistr(x, μ=0, σ=√2)
    return @.  1/(σ*√2π)*exp(-0.5((x-μ)/σ)^2)
end # normaldistr

#------------------------------#
# Vector potential generation  #
#-------------------------------------------------------------------------------
"""
    normal3Donlyz((x0, y0, z0), 
                  (xf, yf, zf), 
                  (nx, ny, nz), 
                  (μx, μy), 
                  (σx, σy),
                  amplitude
                  )
According to given spatial domain, resolution, expectationvalue, standard
deviation and amplitude, creates a vector field who's z-component is normally
distributed in x and y according to the formula fz(i, j) = amplitude *
fx(i)fy(j), where fx and fy are normal distributions in x, and y with
expectation value and std equal to μx, μy, σx, σy, respectively. 
"""
function normal3Donlyz((x0, y0, z0)::Tuple{wpFloat, wpFloat, wpFloat},
                       (xf, yf, zf)::Tuple{wpFloat, wpFloat, wpFloat},
                       (nx, ny, nz)::Tuple{wpInt, wpInt, wpInt},
                       (μx, μy)    ::Tuple{wpFloat, wpFloat},
                       (σx, σy)    ::Tuple{wpFloat, wpFloat},
                       amplitude   ::wpFloat
                       )
    # Create spatial axes and find the grid sizes
    xx = collect(LinRange(x0, xf, nx))
    dx = xx[2] - xx[1]
    yy = collect(LinRange(y0, yf, ny))
    dy = yy[2] - yy[1]
    # Account for single point in the z-axis
    if nz == 1
        zz = [z0]
        dz = 0.
            else
        zz = collect(LinRange(z0, zf, nz))
        dz = zz[2] - zz[1]
    end
    # Initialise the vector field
    ndims = 3
    A = zeros(ndims, nx, ny, nz)
    # Evaluate the z-component of the vecor field to be normally distributed in
    # the x and y dimensions.
    for i = 1:nx
        for j = 1:ny
            A[3,i,j,:] .= amplitude * 
                normaldistr(xx[i], μx, σx) * normaldistr(yy[j], μy, σy)
        end
    end
    return (xx, yy, zz), (dx, dy, dz), A
end # function normal3donlyz


end # module utilities
