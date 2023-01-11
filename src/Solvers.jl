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

# Standard libraries
using LinearAlgebra:    ×
# Internal modules
using WorkingPrecision: wpInt, wpFloat
using Meshes:           Mesh          
using Particles:        specieTable
using Interpolations:   grid
using Schemes:          positionHalfStep

"""
    fullOrbit(pos, vel, specie, bField, eField, dt, scheme)

Solves the Lorentz equation of motion using an arbitrary numerical scheme
(defined by the argument `scheme`).
"""
function fullOrbit(pos         ::Vector{wpFloat},
                   v           ::Vector{wpFloat}, # velocity
                   specie      ::wpInt,
                   mesh        ::Mesh,
                   dt          ::wpFloat,
                   interpolator::Function,
                   scheme      ::Function
                   )
    # Extract particle mass and charge
    m = specieTable[specie, 1]
    q = specieTable[specie, 2]
    # Interpolate fields
    fields, _ = grid(mesh,
                  interpolator,
                  pos)
    B = fields[1]
    E = fields[2]
    #
    a = q/m * (E + v × B) # Accelration by the Lorentz force
    newPos, newVel = scheme(pos, v, a, dt)
    return newPos, newVel
end # funcion fullOrbit


function relFullOrbitExplLeapFrog(pos         ::Vector{wpFloat},
                                  vel         ::Vector{wpFloat}, 
                                  specie      ::wpInt,
                                  mesh        ::Mesh,
                                  dt          ::wpFloat,
                                  interpolator::Function,
                                  scheme      ::Function
                                  )
    # Extract particle mass and charge
    mass  = specieTable[specie, 1]
    charge = specieTable[specie, 2]

    #
    # Step 1: Evaluate half-step in time for position
    #
    posHalf = positionHalfStep(pos, vel, dt)
    # Interpolate fields to this location
    fields, _ = grid(mesh,
                  interpolator,
                  posHalf)
    bField = fields[1]
    eField = fields[2]

    # 
    # Step 2: Evaluate full time step in velocity, which is shceme-dependent.
    #
    velNext = scheme(vel, bField, eField, mass, charge, dt)
    
    #
    # Step 3: Evaluate second half of time step in position
    # 
    posNext = positionHalfStep(posHalf, velNext, dt)

    return posNext, velNext
end # function relFullOrbitExpLeapFrog


function GCA(pos         ::Vector{wpFloat},
             vel         ::Vector{wpFloat},
             specie      ::wpInt,
             mesh        ::Mesh,
             dt          ::wpFloat,
             interpolator::Function,
             scheme      ::Function
             )
    # Interpolate fields to this location
    fields, _ = grid(mesh, interpolator, pos)
    # i, j k are cell corner indexes in mesh. Corresponding to the 
    #  position of the particle.
    bField = fields[1] # The magnetic field
    eField = fields[2] # The electric field
    ∇B = fields[3] # The gradient of the magnetic field.
    

    # Extract particle mass and charge
    m = specieTable[specie, 1]
    q = specieTable[specie, 2]
    
    vperp  = vel[4]    # Particle velocity perpendicular to the magne. field
    vparal = vel[5]    # Particle velocity parallel to the magnetic field
    μ  = vel[6]      # The magnetic moment
    B = norm(bField) # The magnetic field strength
    b̂ = bField/B     # An unit vector pointing in the direction of the
                       #  magnet field
    
    # Electric field component parallel to the magnetic field
    Eparal = eField ⋅ b̂ 
    
    # Compute the acceleration 
    accparal = (q*Eparal - μ*b̂⋅∇B)/m # along the magnetic field lines
    acc = accparal*b̂              # The vector
    # With spatially changing fields, the velocity at this point will not be the
    # same as the last, independent of time.
    velHere = vparal*b̂ + b̂/B × (-c*eField + μ*c/q*∇B)
    
    # Use integration scheme to find velocities at the next time step
    vparalnext    = scheme(vparal, accparal, dt)
    posNext, v = scheme(pos, velHere, acc, dt) # Use v∥next? Will
    # essentially be used if the scheme is euler cromer since a∥ is in acc.
    # norm(velNext) - norm(velhere) should equal v∥next
    #or just? posNext, v = scheme(pos, vel[1:3], acc, dt)
    
    # Compute some auxiliary quantities
    vperpnext = √(norm(v)^2 - vparalnext²) 
    μnext = m*vparalnext²/(2B) #  (should be constant for all times)
    # Maybe μ should be forced constant and kept as a parameter to this solver.
    #   This would require a change in the implementation of solvers and
    #   Patch.run!, where e.g. the particle type is passed to solver. Or that
    #   run! is passed with the particle type, not the patch, such that one may
    #   define different run methods depending on the particle type. 

    velNext = [v[1], v[2], v[3], vperpnext, vparalnext, μnext]
    return posNex, velNext
end # function GCA


function magneticFieldStrengthGradient(mesh, cellIdxI, cellIdxJ, cellIdxK)
    return [0.0, 0.0, 0.0]
end # function magneticFieldStrengthGradient

end # module Solvers
