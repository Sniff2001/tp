#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 Solvers.jl
#
#-------------------------------------------------------------------------------
# Module containing the various solvers.
#-------------------------------------------------------------------------------
module Solvers

# Standard library
using LinearAlgebra: ×, ⋅, norm
# Internal modules
using Particles: specieTable
using Constants: c, cSqrdInv


function fullOrbit(pos, vel, specie, bField, eField, dt, scheme)
    mass   = specieTable[specie, 1]
    charge = specieTable[specie, 2]
    acc = charge/mass * (eField + vel × bField)
    newPos, newVel = scheme(pos, vel, acc, dt)
    return newPos, newVel
end # funcion fullOrbit


"""
    vay(pos, vel, specie, bField, eField, dt, scheme)

Implementation of the relativistic Vay pusher. Adapted from J.-L. Vay (2008)
"""
function vay(pos, vel, specie, bField, eField, dt, scheme)
    mass   = specieTable[specie, 1]
    charge = specieTable[specie, 2]
    # Some factors which use is repeated
    factor1 = charge*dt/2mass
    factor2 = 2factor1

    v = norm(vel)                  # Speed of the particle
    γ = √(1/(1 - v^2*cSqrdInv)) # The relativistic gamma-factor
    u = γ.*vel                     # The relativistic velocity

    # Compute half-step in velocity
    uHalf = u + factor1*(eField + vel × bField)

    # Compute various quantities used to get the full step in velcity
    # See Vay 2008 for details
    uPrime = uHalf + (factor1*eField)
    τ = factor1*bField
    uPrimeNorm = norm(uPrime)
    τNorm = norm(τ)
    uStar = uPrime ⋅ τ/c
    γPrime = √(1 - uPrimeNorm^2*cSqrdInv)
    σ = γPrime^2 - τNorm^2
    # Compute the gamma-factor for full step
    γNext = √(0.5(σ + √(σ^2 + 4(τNorm^2*uPrimeNorm^2))))
    t = τ/γNext
    s = 1/(1 + norm(t)^2)
    # Finally, compute the full step in velocity
    uNext = s*(uPrime + (uPrime ⋅ t)t + uPrime × t)

    # Get the non-relativistic velocity and position
    velNext = uNext/γNext
    posIncrement = velNext*dt  # Simple integration of velocity

    return pos + posIncrement, velNext
    
end # Function vay


end # module Solvers
