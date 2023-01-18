#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 Schemes.jl
#
#-------------------------------------------------------------------------------
# Module containing the various schemes for solving differential equations
#-------------------------------------------------------------------------------


module Schemes

# Standard libraries
using LinearAlgebra:    norm, ×, ⋅
# Internal libraries
using WorkingPrecision: wpFloat, wpInt
using Constants:        c, cSqrdInv

#---------------------------------------#
# Integration of differential equations #
#---------------------------------------#---------------------------------------
function euler(pos::Vector{wpFloat},
               vel::Vector{wpFloat}, 
               acc::Vector{wpFloat},
               dt ::wpFloat
               )
    nextPos = @. pos + vel * dt
    nextVel = @. vel + acc * dt
    return nextPos, nextVel
end # function euler
function euler(vel::wpFloat,
               acc::wpFloat,
               dt ::wpFloat,
               )
    nextVel = @. vel + acc * dt
    return nextVel
end # function euler


function eulerCromer(pos::Vector{wpFloat},
                     vel::Vector{wpFloat}, 
                     acc::Vector{wpFloat},
                     dt ::wpFloat
                     )
    nextVel = @. vel + acc * dt
    nextPos = @. pos + nextVel * dt
    return nextPos, nextVel
end # function eulerCromer
function eulerCromer(vel::wpFloat,
                     acc::wpFloat,
                     dt ::wpFloat
                     )
    nextVel = @. vel + acc * dt
    return nextVel
end # function eulerCromer


function positionHalfStep(pos::Vector{wpFloat},
                          vel::Vector{wpFloat},
                          dt ::wpFloat
                          )
    return pos + 0.5dt*vel
end # function positionHalfStep


"""
    vay(pos, vel, specie, bField, eField, dt, scheme)

Implementation of the Vay pusher. For integrating the relativistic Lorentz
eqution. Adapted from J.-L. Vay (2008). This is only the second step of the
implementation, where the relativistic velocity is advanced a time step. The
first step involves advancing the position half a time step using
`positionHalfStep`, and in the last step the position is advanced its final 
half using the same function but with the new relativistic velocity.
"""
function vay(vel   ::Vector{wpFloat},
             bField::Vector{wpFloat},
             eField::Vector{wpFloat},
             mass  ::wpFloat,
             charge::wpFloat,
             dt    ::wpFloat
             )
    # Some factors which use is repeated
    factor1 = 0.5charge*dt/mass

    v = norm(vel)                  # Speed of the particle
    γ = √(1/(1 - v^2*cSqrdInv)) # The relativistic gamma-factor
    u = γ * vel                   # The relativistic velocity

    # Compute half-step in velocity
    uHalf = u + factor1*(eField + vel × bField)

    # Compute auxiliary quantities used to get the full step in velcity
    # See Vay 2008 for details
    uPrime = uHalf + factor1*eField
    τ = factor1*bField
    uPrimeNorm = norm(uPrime)
    τNorm = norm(τ)
    uStar = uPrime ⋅ τ/c
    γPrime = √(1 + uPrimeNorm^2*cSqrdInv)
    σ = γPrime^2 - τNorm^2
    # Compute the gamma-factor for full step
    γNext = √(0.5(σ + √(σ^2 + 4(τNorm^2 + uStar^2))))
    t = τ/γNext
    s = 1/(1 + norm(t)^2)
    # Finally, compute the full step in velocity
    uNext = s*(uPrime + (uPrime ⋅ t)t + uPrime × t)

    # Get the non-relativistic velocity and position
    velNext = uNext/γNext
    
    return velNext
end # function vay


"""
    boris(pos, vel, specie, bField, eField, dt, scheme)

Implementation of the Boris pusher. For integrating the relativistic Lorentz
eqution. Adapted form Ripperda et al. 2018. This is only the second step of the
implementation, where the relativistic velocity is advanced a time step. The
first step involves advancing the position half a time step using
`positionHalfStep`, and in the last step the position is advanced its final 
half using the same function but with the new relativistic velocity.
"""
function boris(vel   ::Vector{wpFloat},
               bField::Vector{wpFloat},
               eField::Vector{wpFloat},
               mass  ::wpFloat,
               charge::wpFloat,
               dt    ::wpFloat
               )
    # Some factors which use is repeated
    factor1 = 0.5charge*dt/mass

    v = norm(vel)                  # Speed of the particle
    γ = √(1/(1 - v^2*cSqrdInv)) # The relativistic gamma-factor
    u = γ * vel                   # The relativistic velocity
    

    # First half of the electric field acceleration
    uMinus = u + factor1*eField
    uMinusNorm = norm(uMinus)
    
    # Compute auxiliary quantities
    γMinus = √(1 + uMinusNorm*cSqrdInv)
    t = factor1*bField/γMinus
    tNorm = norm(t)
    s = 2t/(1 + tNorm^2)

    # Rotation step
    uPlus = uMinus + (uMinus + (uMinus × t)) × s

    # Second half of electric field acceleration
    uNext = uPlus + factor1*eField
    uNextNorm = norm(uNext)
    γNext = √(1 + uNextNorm^2*cSqrdInv)

    # Get the non-relativistic velocity and position
    velNext = uNext/γNext
    
    return velNext
end # function boris
#-------------------------------------------------------------------------------


#-----------------#
# Differentiation #
#-----------------#-------------------------------------------------------------
"""
    derivateCentral(field, dx)
First and last grid point are ill calculated.
"""
function derivateCentral(field::Vector{wpFloat},
                         dx
                         )
    return (circshift(field, -1) - circshift(field, 1))/2dx
end # function derivateCentral
#|
function derivateCentral(field      ::Array{wpFloat, 3},
                         gridSpacing::wpFloat,
                         axis       ::Tuple{wpInt, wpInt, wpInt}
                         )
    ax1 = -1 .* axis
    return (circshift(field, ax1) - circshift(field, axis))/2gridSpacing
end # function derivateCentral


"""
    derivate4thOrder(field, dx)
First and last two grid point are ill calculated.
"""
function derivate4thOrder(field      ::Array{wpFloat, 3},
                          gridSpacing::wpFloat,
                          axis       ::Tuple{wpInt, wpInt, wpInt}
                          )
    ax1 = -2 .* axis
    ax2 = -1 .* axis
    ax3 =  2 .* axis
    return (-circshift(field, ax1) + 8.0 .* circshift(field, ax2) - 
        8.0 .* circshift(field, axis) + circshift(field, ax3))/12.0gridSpacing
end # function derivate4thOrder


"""
    ∇(field, dx, dy, dz, scheme)
The gradient operator. Calucaletes the gradient of a 3-dimensional scalar
field. Requires grid spacing on all three axis and the numerical scheme as
arguments. The scheme is given as a function type, e.g. Schemes.derivateCentral.
"""
function ∇(field::Array{wpFloat, 3},
           dx   ::wpFloat,
           dy   ::wpFloat,
           dz   ::wpFloat,
           scheme::Function
           )
    ∂f∂x = scheme(field, dx, (1,0,0))
    ∂f∂y = scheme(field, dy, (0,1,0))
    ∂f∂z = scheme(field, dz, (0,0,1))
    nx, ny, nz = size(field)
    gradient = zeros(3, nx, ny, nz)
    gradient[1, :, :, :] = ∂f∂x
    gradient[2, :, :, :] = ∂f∂y
    gradient[3, :, :, :] = ∂f∂z
    return gradient
end # functin ∇ 
#|
function ∇(field::Array{wpFloat, 2},
           dx   ::wpFloat,
           dy   ::wpFloat,
           scheme::Function
           )
    dfdx = scheme(field, dx, (1,0,0))
    dfdy = scheme(field, dy, (0,1,0))
    nx, ny = size(field)
    gradient = zeros(3, nx, ny)
    gradient[1, :, :] = ∂f∂x
    gradient[2, :, :] = ∂f∂y
    return gradient
end # function ∇ 


"""
    ∇×(field, dx, dy, dz, derivscheme)
The curl operator. Calculates the curl of a 3-dimensional vector
field. Requires grid spacing on all three axis and the numerical scheme as
arguments. The scheme is given as a function type, e.g. Schemes.derivateCentral.
"""
function curl(field      ::Array{wpFloat, 4},
              gridsizes  ::Tuple{wpFloat, wpFloat, wpFloat},
              derivscheme::Function
              )
    dx, dy, dz = gridsizes
    fx = field[1,:,:,:]
    fy = field[2,:,:,:]
    fz = field[3,:,:,:]
    derivx = (1,0,0)
    derivy = (0,1,0)
    derivz = (0,0,1)
    ∂fx∂y = derivscheme(fx, dy, derivy)
    ∂fx∂z = derivscheme(fx, dz, derivz)
    ∂fy∂x = derivscheme(fy, dx, derivx)
    ∂fy∂z = derivscheme(fy, dz, derivz)
    ∂fz∂x = derivscheme(fz, dx, derivx)
    ∂fz∂y = derivscheme(fz, dy, derivy)
    _, nx, ny, nz = size(field)
    result = zeros(3, nx, ny, nz)
    result[1, :, :, :] = ∂fz∂y .- ∂fy∂z
    result[2, :, :, :] = ∂fx∂z .- ∂fz∂x
    result[3, :, :, :] = ∂fy∂x .- ∂fx∂y
    return result
end # functin curl


end # module schemes
