#-------------------------------------------------------------------------------
# Created 05.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                Particles.jl
#
#-------------------------------------------------------------------------------
# Module containing the Particle structs and methods
#-------------------------------------------------------------------------------

module Particles

using WorkingPrecision: wpFloat, wpInt
using Constants: m_e, m_p, e
using Utilities: norm3


#-------------#   
# Exports     # 
#-------------#-----------------------------------------------------------------
export TraceParticle
export ParticleSoA # Particles represented as struct of arrays
export GCAParticleSoA
export specieTable # Maping specie to mass and charge
export getpos
export getvel
export reset!      # Resets particle positions to zero (except initial position)
export setinitpos! # Sets the initial position of particles
export setinitvel! # Sets the initial velocity of particles
export kineticenergy # Computes the non-rel. kinetic energy at all time steps

#------------------#
# Global variables #
#------------------#------------------------------------------------------------
#              mass charge
specieTable = [m_e  -e     # Electron
               m_p   e     # Proton
               1.0 3.0     # Unit mass and charge = 3
               1.0 1.0]    # Unit mass and charge

#-------------#   
# Structs     # 
#-------------#-----------------------------------------------------------------
"""
    TraceParticle
The supertype of all trace particles
"""
abstract type TraceParticle end

    
mutable struct ParticleSoA <: TraceParticle
    pos    ::Array{wpFloat, 3}
    vel    ::Array{wpFloat, 3}
    species::Vector{wpInt}   # Particle specie identifier (e.g. electron, proton)
    alive  ::Vector{Bool}
    weight ::Vector{wpFloat}
    
    
    # Constructors
    #--------------------------------------------------------------------------
    """
        ParticleSoA(pos::Matrix, vel::Matrix, specie, numSteps)

    One would normally only pass initial conditions. This constructor handles
    the creation of the type accordingly, by adding the initial conditions to
    an higher order array.
    """
    function ParticleSoA(pos     ::Matrix{wpFloat},
                         vel     ::Matrix{wpFloat},
                         species ::Vector{wpInt},
                         numSteps::Integer
                         )
        numDims, numParticles = size(pos)
        numVels, numParticles = size(vel)
        positions  = zeros(wpFloat, numDims, numParticles, numSteps + 1)
        velocities = zeros(wpFloat, numVels, numParticles, numSteps + 1)
        positions[:, :, 1] .= pos
        velocities[:, :, 1] .= vel
        alive = ones(Bool, numParticles)
        weight = ones(wpFloat, numParticles)
        return new(positions, velocities, species, alive, weight)
    end # constructor 
end # mutable struct ParticleSoA


mutable struct GCAParticleSoA <: TraceParticle
    R      ::Array{wpFloat, 3} # The position of the guiding centre
    vparal ::Matrix{wpFloat}   # The velocity parallel to the magnetic field
    μ      ::Vector{wpFloat}   # Magnetic moment μ of particle
    species::Vector{wpInt}     # Particle specie identifier (e.g. electron, proton)
    alive  ::Vector{Bool}
    weight ::Vector{wpFloat}
    
    
    # Constructors
    #--------------------------------------------------------------------------
    """
        GCAParticleSoA(R, vparal, μ, species, alive, weight)

    One would normally only pass initial conditions. This constructor handles
    the creation of the type accordingly, by adding the initial conditions to
    an higher order array.
    """
    function GCAParticleSoA(
        initR     ::Matrix{wpFloat},
        initvparal::Vector{wpFloat},
        μ         ::Vector{wpFloat},
        species   ::Vector{wpInt},
        numSteps  ::Integer
        )
        numDims, numParticles = size(initR)
        R      = zeros(wpFloat, numDims, numParticles, numSteps + 1)
        vparal = zeros(wpFloat, numParticles, numSteps + 1)
        R[:, :, 1] .= initR
        vparal[:, :, 1] .= initvparal
        alive = ones(Bool, numParticles)
        weight = ones(wpFloat, numParticles)
        return new(R, vparal, μ, species, alive, weight)
    end # constructor 
end # mutable struct ParticleSoA

#------------------------#
# Particle get-functions #
#-------------------------------------------------------------------------------
function getpos(particles::ParticleSoA)
    return particles.pos
end
#|
function getpos(particles::GCAParticleSoA)
    return particles.R
end


function getvel(particles::ParticleSoA)
    return particles.vel
end
#------------------------#
# Particle set-functions #
#-------------------------------------------------------------------------------
function reset!(particles::ParticleSoA)
    n = length(particles.pos[1,1,:])
    particles.pos[:, :, 2:n] .= 0.0
    particles.vel[:, :, 2:n] .= 0.0
end #function reset!


function setinitpos!(particles::ParticleSoA,
                     pos      ::Matrix{wpFloat})
    particles.pos[:, :, 1] .= pos
end # function setinitpos
#|
function setinitpos!(particles::ParticleSoA,
                     pos      ::Vector{wpFloat},
                     partIdx  ::wpInt)
    particles.pos[:, partIdx, 1] .= pos
end # function setinitpos


function setinitvel!(particles::ParticleSoA,
                     vel      ::Matrix{wpFloat})
    particles.vel[:, :, 1] .= vel
end # function setinitvel
#|
function setinitvel!(particles::ParticleSoA,
                     vel      ::Vector{wpFloat},
                     partIdx  ::wpInt)
    particles.vel[:, partIdx, 1] .= vel
end # function setinitvel


#----------------------#
# Auxiliary quantities #
#-------------------------------------------------------------------------------
function kineticenergy(particles::ParticleSoA)
    _, npart, N = size(particles.pos) 
    v = norm3(particles.vel) 
    Ek = zeros(npart, N)
    for i = 1:npart
        mass = specieTable[particles.species[i], 1]
        @. Ek[i,:] = 0.5*mass*v[i,:]^2
    end # loop i
    return Ek
end # function kineticenergy
#|
function kineticenergy(particles::GCAParticleSoA)
    npart, N = size(particles.vparal) 
    v = particles.vparal
    Ek = zeros(npart, N)
    for i = 1:npart
        mass = specieTable[particles.species[i], 1]
        @. Ek[i,:] = 0.5*mass*v[i,:]^2
    end # loop i
    return Ek
end # function kineticenergy

end # module particles
