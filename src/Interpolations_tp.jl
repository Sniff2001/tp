#-------------------------------------------------------------------------------
# Created 12.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#             Interpolations_tp.jl
#
#-------------------------------------------------------------------------------
# Module containing interpolation methods.
#-------------------------------------------------------------------------------

module Interpolations_tp

using Meshes
using Interpolations

export gridinterp
export locateCell
export trilinear
export bilinear_xz


#-----------------------#
# Type of interpolation #
#-------------------------------------------------------------------------------
function gridinterp(
    tensorfield ::AbstractArray{T} where {T<:Real},
    interpolator::Function,
    pos         ::Vector{T} where {T<:Real},
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real}
    )
    numDims = 3
    x, y, z = pos
    # Find which cell the position is in, and store the indices of its corner.
    i = locateCell(xx, x)
    j = locateCell(yy, y)
    k = locateCell(zz, z)
    interpolatedfield = interpolator(
        tensorfield,
        xx, yy, zz,
        (i,j,k),
        (x,y,z)
    )
    return interpolatedfield, i, j, k
end # function gridinterp
#|
function gridinterp(
    mesh        ::Mesh,
    interpolator::Function,
    pos         ::Vector{T} where {T<:Real}
    )
    numDims = 3
    x, y, z = pos
    # Find which cell the position is in, and store the indices of its corner.
    i = locateCell(mesh.xCoords, x)
    j = locateCell(mesh.yCoords, y)
    k = locateCell(mesh.zCoords, z)
    # Interpolate
    interpolatedfields = interpolator(mesh, (i,j,k), (x,y,z))
    return interpolatedfields, i, j, k
end # function gridinterp



#-----------#
# Utilities #
#-------------------------------------------------------------------------------
"""
    locateCell(coords, point)

Finds the position of `point` in the vector `coords` using binary search.
Returns its lower neighbour.
"""
function locateCell(coords::Vector{T} where {T<:Real},
                    point ::Number
                    )
    # Initial bounds
    low = 1
    high = length(coords)
    while  high > low + 1 # While we haven't found the cell
        mid  = floor(Int64, (high + low)/2)
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


#---------------------#
# Interpolation steps #
#-------------------------------------------------------------------------------
function trilinear( # Arbitrary vector-field
    tensorfield ::AbstractArray{T} where {T<:Real},
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real},
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{Real, Real, Real}
    )
    coefficients = trilinearcoefficients(xx,
                                         yy,
                                         zz,
                                         (i,j,k),
                                         (x,y,z))
    f = trilinearsum(tensorfield,
                     (i,j,k),
                     coefficients)
    return f
end # function trilinear
#|
function trilinear( # Method for passing the mesh-struct
    mesh   ::Mesh,
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{Real, Real, Real}
    )
    coefficients = trilinearcoefficients(mesh.xCoords,
                                         mesh.yCoords,
                                         mesh.zCoords,
                                         (i,j,k),
                                         (x,y,z))
    B = trilinearsum(mesh.bField,
                              (i,j,k),
                              coefficients)
    E = trilinearsum(mesh.eField,
                              (i,j,k),
                              coefficients)
    return B, E
end # function trilinear


function bilinear_xz( # Arbitrary vector-field
    tensorfield ::AbstractArray{T} where {T<:Real},
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real},
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{Real, Real, Real}
    )
    coefficients = bilinearcoefficients(xx,
                                        zz,
                                        (i,j),
                                        (x,z))
    f = bilinearsum(dropdims(tensorfield, dims=ndims(tensorfield)-1),
                    (i,j),
                    coefficients)
    return f
end # function trilinear
#|
function bilinear_xz( 
    mesh   ::Mesh,
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{Real, Real, Real}
    )
    coefficients = bilinearcoefficients(mesh.xCoords,
                                        mesh.zCoords,
                                        (i,j),
                                        (x,z))
    B = bilinearsum(dropdims(mesh.bField, dims=3),
                    (i,j),
                    coefficients)
    E = bilinearsum(dropdims(mesh.eField, dims=3),
                    (i,j),
                    coefficients)
    return B, E
end # function trilinear


function trilinearGCA(
    mesh   ::Mesh,
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{Real, Real, Real}
    )
    coefficients = trilinearcoefficients(mesh.xCoords,
                                         mesh.yCoords,
                                         mesh.xCoords,
                                         (i,j,k),
                                         (x,y,z))
    B = trilinearsum(mesh.bField,
                              (i,j,k),
                              coefficients)
    E = trilinearsum(mesh.eField,
                              (i,j,k),
                              coefficients)
    ∇B = trilinearsum(mesh.∇B,
                              (i,j,k),
                              coefficients)
    ∇b̂ = trilinearsum(mesh.∇b̂,
                              (i,j,k),
                              coefficients)
    ∇ExB = trilinearsum(mesh.∇ExB,
                              (i,j,k),
                              coefficients)

    return B, E, ∇B, ∇b̂, ∇ExB
end # function trilinearGCA


function trilinearcoefficients(
    xCoords::Vector{T} where {T<:Real},
    yCoords::Vector{T} where {T<:Real},
    zCoords::Vector{T} where {T<:Real},
    (i,j,k)::Tuple{Integer, Integer, Integer},
    (x,y,z)::Tuple{Real, Real, Real}
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
end # function trinlinearcoefficients


function bilinearcoefficients(
    xCoords::Vector{T} where {T<:Real},
    yCoords::Vector{T} where {T<:Real},
    (i,j)::Tuple{Integer, Integer},
    (x,y)::Tuple{Real, Real}
    )
    t = (x - xCoords[i])/(xCoords[i + 1] - xCoords[i])
    u = (y - yCoords[j])/(yCoords[j + 1] - yCoords[j])
    #
    c0 = (1 - t)*(1 - u)
    c1 = t*(1 - u)
    c2 = t*u
    c3 = (1 - t)*u
    #
    return c0, c1, c2, c3
end # bilinearcoefficients

function trilinearsum(
    tensorfield::Array{T, 5} where {T<:Real},
    (i,j,k)    ::Tuple{Integer, Integer, Integer},
    c          ::NTuple{8, Real}
    )
    c0, c1, c2, c3, c4, c5, c6, c7 = c
    A0 = tensorfield[:, :,   i,   j,   k]
    A1 = tensorfield[:, :, i+1,   j,   k]
    A2 = tensorfield[:, :, i+1, j+1,   k]
    A3 = tensorfield[:, :,   i, j+1,   k]
    A4 = tensorfield[:, :,   i,   j, k+1]
    A5 = tensorfield[:, :, i+1,   j, k+1]
    A6 = tensorfield[:, :, i+1, j+1, k+1]
    A7 = tensorfield[:, :,   i, j+1, k+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 + c4*A4 + c5*A5 + c6*A6 + c7*A7
    return A
end # function trilinearsum
#|
function trilinearsum(
    vectorfield::Array{T, 4} where {T<:Real},
    (i,j,k)    ::Tuple{Integer, Integer, Integer},
    c          ::NTuple{8, Real}
    )
    c0, c1, c2, c3, c4, c5, c6, c7 = c
    A0 = vectorfield[:,   i,   j,   k]
    A1 = vectorfield[:, i+1,   j,   k]
    A2 = vectorfield[:, i+1, j+1,   k]
    A3 = vectorfield[:,   i, j+1,   k]
    A4 = vectorfield[:,   i,   j, k+1]
    A5 = vectorfield[:, i+1,   j, k+1]
    A6 = vectorfield[:, i+1, j+1, k+1]
    A7 = vectorfield[:,   i, j+1, k+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 + c4*A4 + c5*A5 + c6*A6 + c7*A7
    return A
end # function trilinearsum
#|
function trilinearsum(
    scalarfield::Array{T, 3} where {T<:Real},
    (i,j,k)    ::Tuple{Integer, Integer, Integer},
    c          ::NTuple{8, Real}
    )
    c0, c1, c2, c3, c4, c5, c6, c7 = c
    A0 = scalarfield[  i,   j,   k]
    A1 = scalarfield[i+1,   j,   k]
    A2 = scalarfield[i+1, j+1,   k]
    A3 = scalarfield[  i, j+1,   k]
    A4 = scalarfield[  i,   j, k+1]
    A5 = scalarfield[i+1,   j, k+1]
    A6 = scalarfield[i+1, j+1, k+1]
    A7 = scalarfield[  i, j+1, k+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 + c4*A4 + c5*A5 + c6*A6 + c7*A7
    return A
end # function trilinearsum


function bilinearsum(
    tensorfield::Array{T, 4} where {T<:Real},
    (i,j)      ::Tuple{Integer, Integer},
    c          ::NTuple{4, Real}
    )
    c0, c1, c2, c3 = c
    A0 = tensorfield[:, :,   i,   j]
    A1 = tensorfield[:, :, i+1,   j]
    A2 = tensorfield[:, :, i+1, j+1]
    A3 = tensorfield[:, :,   i, j+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 
    return A
end # bilinearsum
#|
function bilinearsum(
    vectorfield::Array{T, 3} where {T<:Real},
    (i,j)      ::Tuple{Integer, Integer},
    c          ::NTuple{4, Real}
    )
    c0, c1, c2, c3 = c
    A0 = vectorfield[:,   i,   j]
    A1 = vectorfield[:, i+1,   j]
    A2 = vectorfield[:, i+1, j+1]
    A3 = vectorfield[:,   i, j+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 
    return A
end # bilinearsum
#|
function bilinearsum(
    scalarfield::Array{T, 2} where {T<:Real},
    (i,j)      ::Tuple{Integer, Integer},
    c          ::NTuple{4, Real}
    )
    c0, c1, c2, c3 = c
    A0 = scalarfield[  i,   j]
    A1 = scalarfield[i+1,   j]
    A2 = scalarfield[i+1, j+1]
    A3 = scalarfield[  i, j+1]
    A  = c0*A0 + c1*A1 + c2*A2 + c3*A3 
    return A
end # bilinearsum

function trilinear_ip(
    tensorfield1 ::AbstractArray{T} where {T<:Real},
    tensorfield2 ::AbstractArray{T} where {T<:Real},
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real},
    pbc         ::Tuple{Bool, Bool, Bool}
    )
    tensorfield1_interp_x = interpolate((xx,yy,zz), tensorfield1[1,:,:,:], Gridded(Linear()))
    tensorfield1_interp_y = interpolate((xx,yy,zz), tensorfield1[2,:,:,:], Gridded(Linear()))
    tensorfield1_interp_z = interpolate((xx,yy,zz), tensorfield1[3,:,:,:], Gridded(Linear()))
    tensorfield2_interp_x = interpolate((xx,yy,zz), tensorfield2[1,:,:,:], Gridded(Linear()))
    tensorfield2_interp_y = interpolate((xx,yy,zz), tensorfield2[2,:,:,:], Gridded(Linear()))
    tensorfield2_interp_z = interpolate((xx,yy,zz), tensorfield2[3,:,:,:], Gridded(Linear()))
    
    if pbc[1]
        tensorfield1_interp_x = extrapolate(tensorfield1_interp_x, Periodic())
        tensorfield2_interp_x = extrapolate(tensorfield2_interp_x, Periodic())
    else
        tensorfield1_interp_x = extrapolate(tensorfield1_interp_x, Flat())
        tensorfield2_interp_x = extrapolate(tensorfield2_interp_x, Flat())
    end

    if pbc[2]
        tensorfield1_interp_y = extrapolate(tensorfield1_interp_y, Periodic())
        tensorfield2_interp_y = extrapolate(tensorfield2_interp_y, Periodic())
    else
        tensorfield1_interp_y = extrapolate(tensorfield1_interp_y, Flat())
        tensorfield2_interp_y = extrapolate(tensorfield2_interp_y, Flat())
    end

    if pbc[3]
        tensorfield1_interp_z = extrapolate(tensorfield1_interp_z, Periodic())
        tensorfield2_interp_z = extrapolate(tensorfield2_interp_z, Periodic())
    else
        tensorfield1_interp_z = extrapolate(tensorfield1_interp_z, Flat())
        tensorfield2_interp_z = extrapolate(tensorfield2_interp_z, Flat())
    end

    tensorfield1_interp = (x, y, z) -> [tensorfield1_interp_x(x,y,z), tensorfield1_interp_y(x,y,z), tensorfield1_interp_z(x,y,z)]
    tensorfield2_interp = (x, y, z) -> [tensorfield2_interp_x(x,y,z), tensorfield2_interp_y(x,y,z), tensorfield2_interp_z(x,y,z)]
    return (tensorfield1_interp, tensorfield2_interp)
end # trilinear_ip

# currently i doubt bsplines work, WIP
function quadratic_bspline(
    tensorfield1 ::AbstractArray{T} where {T<:Real},
    tensorfield2 ::AbstractArray{T} where {T<:Real},
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real},
    pbc         ::Tuple{Bool, Bool, Bool}
    )
    tensorfield1_interp_x = interpolate(tensorfield1[1,:,:,:], BSpline(Quadratic()))
    tensorfield1_interp_y = interpolate(tensorfield1[2,:,:,:], BSpline(Quadratic()))
    tensorfield1_interp_z = interpolate(tensorfield1[3,:,:,:], BSpline(Quadratic()))
    tensorfield2_interp_x = interpolate(tensorfield2[1,:,:,:], BSpline(Quadratic()))
    tensorfield2_interp_y = interpolate(tensorfield2[2,:,:,:], BSpline(Quadratic()))
    tensorfield2_interp_z = interpolate(tensorfield2[3,:,:,:], BSpline(Quadratic()))

    xa = LinRange(xx[1], xx[end], length(xx))
    ya = LinRange(yy[1], yy[end], length(yy))
    za = LinRange(zz[1], zz[end], length(zz))

    tensorfield1_scaledinterp_x = scale(tensorfield1_interp_x, xa, ya, za)
    tensorfield2_scaledinterp_x = scale(tensorfield2_interp_x, xa, ya, za)
    tensorfield1_scaledinterp_y = scale(tensorfield1_interp_y, xa, ya, za)
    tensorfield2_scaledinterp_y = scale(tensorfield2_interp_y, xa, ya, za)
    tensorfield1_scaledinterp_z = scale(tensorfield1_interp_z, xa, ya, za)
    tensorfield2_scaledinterp_z = scale(tensorfield2_interp_z, xa, ya, za)

    if pbc[1]
        tensorfield1_scaledinterp_x = extrapolate(tensorfield1_scaledinterp_x, Periodic())
        tensorfield2_scaledinterp_x = extrapolate(tensorfield2_scaledinterp_x, Periodic())
    else
        tensorfield1_scaledinterp_x = extrapolate(tensorfield1_scaledinterp_x, Flat())
        tensorfield2_scaledinterp_x = extrapolate(tensorfield2_scaledinterp_x, Flat())
    end

    if pbc[2]
        tensorfield1_scaledinterp_y = extrapolate(tensorfield1_scaledinterp_y, Periodic())
        tensorfield2_scaledinterp_y = extrapolate(tensorfield2_scaledinterp_y, Periodic())
    else
        tensorfield1_scaledinterp_y = extrapolate(tensorfield1_scaledinterp_y, Flat())
        tensorfield2_scaledinterp_y = extrapolate(tensorfield2_scaledinterp_y, Flat())
    end

    if pbc[3]
        tensorfield1_scaledinterp_z = extrapolate(tensorfield1_scaledinterp_z, Periodic())
        tensorfield2_scaledinterp_z = extrapolate(tensorfield2_scaledinterp_z, Periodic())
    else
        tensorfield1_scaledinterp_z = extrapolate(tensorfield1_scaledinterp_z, Flat())
        tensorfield2_scaledinterp_z = extrapolate(tensorfield2_scaledinterp_z, Flat())
    end

    tensorfield1_interp = (x, y, z) -> [tensorfield1_scaledinterp_x(x,y,z), tensorfield1_scaledinterp_y(x,y,z), tensorfield1_scaledinterp_z(x,y,z)]
    tensorfield2_interp = (x, y, z) -> [tensorfield2_scaledinterp_x(x,y,z), tensorfield2_scaledinterp_y(x,y,z), tensorfield2_scaledinterp_z(x,y,z)]
    return (tensorfield1_interp, tensorfield2_interp)
end # quadratic_bspline

function cubic_bspline(
    tensorfield1 ::AbstractArray{T} where {T<:Real},
    tensorfield2 ::AbstractArray{T} where {T<:Real},
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real},
    pbc         ::Tuple{Bool, Bool, Bool}
    )
    tensorfield1_interp_x = interpolate(tensorfield1[1,:,:,:], BSpline(Cubic(Line(OnGrid()))))
    tensorfield1_interp_y = interpolate(tensorfield1[2,:,:,:], BSpline(Cubic(Line(OnGrid()))))
    tensorfield1_interp_z = interpolate(tensorfield1[3,:,:,:], BSpline(Cubic(Line(OnGrid()))))
    tensorfield2_interp_x = interpolate(tensorfield2[1,:,:,:], BSpline(Cubic(Line(OnGrid()))))
    tensorfield2_interp_y = interpolate(tensorfield2[2,:,:,:], BSpline(Cubic(Line(OnGrid()))))
    tensorfield2_interp_z = interpolate(tensorfield2[3,:,:,:], BSpline(Cubic(Line(OnGrid()))))

    xa = LinRange(xx[1], xx[end], length(xx))
    ya = LinRange(yy[1], yy[end], length(yy))
    za = LinRange(zz[1], zz[end], length(zz))

    tensorfield1_scaledinterp_x = scale(tensorfield1_interp_x, xa, ya, za)
    tensorfield2_scaledinterp_x = scale(tensorfield2_interp_x, xa, ya, za)
    tensorfield1_scaledinterp_y = scale(tensorfield1_interp_y, xa, ya, za)
    tensorfield2_scaledinterp_y = scale(tensorfield2_interp_y, xa, ya, za)
    tensorfield1_scaledinterp_z = scale(tensorfield1_interp_z, xa, ya, za)
    tensorfield2_scaledinterp_z = scale(tensorfield2_interp_z, xa, ya, za)

    if pbc[1]
        tensorfield1_scaledinterp_x = extrapolate(tensorfield1_scaledinterp_x, Periodic())
        tensorfield2_scaledinterp_x = extrapolate(tensorfield2_scaledinterp_x, Periodic())
    else
        tensorfield1_scaledinterp_x = extrapolate(tensorfield1_scaledinterp_x, Flat())
        tensorfield2_scaledinterp_x = extrapolate(tensorfield2_scaledinterp_x, Flat())
    end

    if pbc[2]
        tensorfield1_scaledinterp_y = extrapolate(tensorfield1_scaledinterp_y, Periodic())
        tensorfield2_scaledinterp_y = extrapolate(tensorfield2_scaledinterp_y, Periodic())
    else
        tensorfield1_scaledinterp_y = extrapolate(tensorfield1_scaledinterp_y, Flat())
        tensorfield2_scaledinterp_y = extrapolate(tensorfield2_scaledinterp_y, Flat())
    end

    if pbc[3]
        tensorfield1_scaledinterp_z = extrapolate(tensorfield1_scaledinterp_z, Periodic())
        tensorfield2_scaledinterp_z = extrapolate(tensorfield2_scaledinterp_z, Periodic())
    else
        tensorfield1_scaledinterp_z = extrapolate(tensorfield1_scaledinterp_z, Flat())
        tensorfield2_scaledinterp_z = extrapolate(tensorfield2_scaledinterp_z, Flat())
    end

    tensorfield1_interp = (x, y, z) -> [tensorfield1_scaledinterp_x(x,y,z), tensorfield1_scaledinterp_y(x,y,z), tensorfield1_scaledinterp_z(x,y,z)]
    tensorfield2_interp = (x, y, z) -> [tensorfield2_scaledinterp_x(x,y,z), tensorfield2_scaledinterp_y(x,y,z), tensorfield2_scaledinterp_z(x,y,z)]
    return (tensorfield1_interp, tensorfield2_interp)
end # quadratic_bspline

end # module Interpolations_tp

