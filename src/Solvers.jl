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
using LinearAlgebra:     ×, ⋅, norm
# Internal modules
using Constants:         c
using Meshes:            Mesh          
using Particles:         specieTable
using Interpolations_tp: gridinterp
using Schemes:           positionHalfStep

"""
    fullOrbit(pos, vel, specie, bField, eField, dt, scheme)

Solves the Lorentz equation of motion using an arbitrary numerical scheme
(defined by the argument `scheme`).
"""
function fullOrbit_interstaticfield(
    pos         ::Vector{T} where {T<:Real},
    v           ::Vector{T} where {T<:Real}, # velocity
    specie      ::Integer,
    mesh        ::Mesh,
    dt          ::Real,
    interpolator::Function,
    scheme      ::Function
    )
    # Extract particle mass and charge
    m = specieTable[specie, 1]
    q = specieTable[specie, 2]
    # Interpolate fields
    fields, _ = gridinterp(mesh,
                           interpolator,
                           pos)
    B = fields[1]
    E = fields[2]
    #
    statevector = [pos; v]
    statevectorNext = scheme(statevector,
                             dt, 
                             eomLorentzforce,
                             B, E, q, m)
    return statevectorNext[1:3], statevectorNext[4:6]
end # funcion fullOrbit_interstaticfield

function fullOrbit(pos        ::Vector{T} where {T<:Real},
                   v           ::Vector{T} where {T<:Real}, # velocity
                   specie      ::Integer,
                   mesh        ::Mesh,
                   dt          ::Real,
                   interpolator::Function,
                   scheme      ::Function
                   )
    # Extract particle mass and charge
    m = specieTable[specie, 1]
    q = specieTable[specie, 2]
    #
    statevector = [pos; v]
    statevectorNext = scheme(statevector,
                             dt, 
                             eomLorentzforce,
                             q, m, mesh, interpolator,
                             )
    return statevectorNext[1:3], statevectorNext[4:6]
end # function fullOrbit

function fullOrbit(
    pos         ::Vector{T} where {T<:Real},
    v           ::Vector{T} where {T<:Real}, # velocity
    specie      ::Integer,
    mesh        ::Mesh,
    dt          ::Real,
    interpolator::Tuple,
    scheme      ::Function
    )
# Extract particle mass and charge
m = specieTable[specie, 1]
q = specieTable[specie, 2]
B = interpolator[1]
E = interpolator[2]
#
statevector = [pos; v]
statevectorNext = scheme(statevector,
              dt, 
              eomLorentzforce,
              B, E, q, m
              )
return statevectorNext[1:3], statevectorNext[4:6]
end # function fullOrbit


function eomLorentzforce(
    s::Vector{T} where {T<:Real}, # The state vector
    B::Vector{T} where {T<:Real}, # The magnetic field
    E::Vector{T} where {T<:Real}, # The electric field
    q::Real,         # Charge
    m::Real          # Mass
    )
    x = s[1:3] # The position vector
    v = s[4:6] # The velocity vector
    dvdt = q/m * (E + v × B) 
    dxdt = v
    dsdt = [dxdt; dvdt]
    return dsdt
end # function eomLorentzforce
#|
function eomLorentzforce(
    s::Vector{T} where {T<:Real}, # The state vector
    B::Function, # The magnetic field
    E::Function, # The electric field
    q::Real,         # Charge
    m::Real          # Mass
    )
    x = s[1:3] # The position vector
    v = s[4:6] # The velocity vector
    dvdt = q/m * (E(x[1], x[2], x[3]) + v × B(x[1], x[2], x[3])) 
    dxdt = v
    dsdt = [dxdt; dvdt]
    return dsdt
end # function eomLorentzforce
#|
function eomLorentzforce(
    statevector ::Vector{T} where {T<:Real}, # The state vector
    q           ::Real,         # Charge
    m           ::Real,         # Mass
    mesh        ::Mesh,            # The mesh containing the magnetic field
    interpolator::Function         # Interpolation function used for evaluating
                                   #   the field at the stavector-location
    )
    x = statevector[1:mesh.numdims] # The position vector
    v = statevector[mesh.numdims+1:2mesh.numdims] # The velocity vector
    # Interpolate fields
    fields, _ = gridinterp(mesh,
                           interpolator,
                           x)
    B = fields[1]
    E = fields[2]
    dvdt = q/m * (E + v × B)
    dxdt = v
    dsdt = [dxdt; dvdt]
    return dsdt
end # function eomLorentzforce


function fieldtracingforward(
    statevector ::Vector{T} where {T<:Real}, # Should be just position
    vectorfield ::Array{T, 4} where {T<:Real},
    interpolator::Function,
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real}
    )
    interpfield, _ = gridinterp(vectorfield, interpolator, statevector,
                                xx, yy, zz)
    fieldstrength = norm(interpfield)
    fielddirection = interpfield ./ fieldstrength
    dsdt = fielddirection
    return dsdt
end # function fieldtracing

function fieldtracingbackward(
    statevector ::Vector{T} where {T<:Real}, # Should be just position
    vectorfield ::Array{T, 4} where {T<:Real},
    interpolator::Function,
    xx          ::Vector{T} where {T<:Real},
    yy          ::Vector{T} where {T<:Real},
    zz          ::Vector{T} where {T<:Real}
    )
    interpfield, _ = gridinterp(vectorfield, interpolator, statevector,
                                xx, yy, zz)
    fieldstrength = norm(interpfield)
    fielddirection = interpfield ./ fieldstrength
    dsdt = -fielddirection
    return dsdt
end # function fieldtracing
    

function relFullOrbitExplLeapFrog(pos         ::Vector{T} where {T<:Real},
                                  vel         ::Vector{T} where {T<:Real}, 
                                  specie      ::Integer,
                                  mesh        ::Mesh,
                                  dt          ::Real,
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
    fields, _ = gridinterp(mesh,
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


function GCA(pos         ::Vector{T} where {T<:Real},
             vel         ::Vector{T} where {T<:Real},
             specie      ::Integer,
             mesh        ::Mesh,
             dt          ::Real,
             interpolator::Function,
             scheme      ::Function
             )
    # Interpolate fields to this location
    fields, _ = gridinterp(mesh, interpolator, pos)
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
    vperpnext = √(norm(v)^2 - vparalnext^2) 
    μnext = m*vperpnext^2/(2B) #  (should be constant for all times)
    # Maybe μ should be forced constant and kept as a parameter to this solver.
    #   This would require a change in the implementation of solvers and
    #   Patch.run!, where e.g. the particle type is passed to solver. Or that
    #   run! is passed with the particle type, not the patch, such that one may
    #   define different run methods depending on the particle type. 

    velNext = [v[1], v[2], v[3], vperpnext, vparalnext, μnext]
    return posNext, velNext
end # function GCA

function GCA(
    pos         ::Vector{T} where {T<:Real},
    vel         ::Real,
    μ           ::Real,
    specie      ::Integer,
    mesh        ::Mesh,
    dt          ::Real,
    interpolator::Function,
    scheme      ::Function
    )
    # Extract particle mass and charge
    m = specieTable[specie, 1]
    q = specieTable[specie, 2]
    #
    statevector = [pos; vel]
    statevectorNext = scheme(statevector,
                             dt, 
                             eomGCA,
                             q, m, μ, mesh, interpolator
                             )
    return statevectorNext[1:3], statevectorNext[4]
    
end # function GCA

function eomGCA(
    statevector ::Vector{T} where {T<:Real},
    q           ::Real,
    m           ::Real,
    μ           ::Real,
    mesh        ::Mesh,
    interpolator::Function
    )
    R      = statevector[1:3]
    vparal = statevector[4] # Particle velocity parallel to the magnetic field

    # Interpolate fields to this location
    fields, _ = gridinterp(mesh, interpolator, R)
    # i, j k are cell corner indexes in mesh. Corresponding to the 
    #  position of the particle.
    B⃗ = fields[1] # The magnetic field
    E⃗ = fields[2] # The electric field
    ∇B = fields[3]     # The gradient of the magnetic field.
    ∇b̂ = fields[4]
    ∇ExB = fields[5]
    local B = norm(B⃗)   # The magnetic field strength
    b̂ = B⃗/B       # An unit vector pointing in the direction of the
                       #  magnetic field
    # Electric field component parallel to the magnetic field
    Eparal = E⃗⋅b̂ 
    # Calculate drifts
    ExBdrift = (E⃗ × b̂)/B
    ∇Bdrift = μ/(q*B)*(b̂ × ∇B)
    # Total time derivatives. Assumes ∂/∂t = 0,
    db̂dt = vparal * (∇b̂ * b̂) + ∇b̂*ExBdrift
    dExBdt = vparal * (∇ExB * b̂) + ∇ExB*ExBdrift
    
    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + m*b̂/(q*B) × (vparal*db̂dt + dExBdt)
    #dRperpdt = b̂/B × (-E⃗ + μ/q*∇B + m/q*(vparal*db̂dt + dExBdt))  # old

    # Compute the acceleration 
    dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m # along the magnetic field lines
    # With correction proposed by Birn et al., 2004:
    #dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m + ExBdrift⋅db̂dt + ∇Bdrift⋅db̂dt
    #dRperpdt = b̂/B × (-c*E⃗ + μ*c/q * ∇B) #old

    # Compute the velocity
    dRdt = vparal*b̂ + dRperpdt
    
    # How to store the perpendicular velocity? Would need to know b̂ at
    #   each R calculate vperp at each R. Could be interesting to store this
    #   as an auxiliary variable somehow, to se how the drift evolves.
    dsdt = [dRdt; dvparaldt]
    return dsdt
end # function GCA


end # module Solvers
