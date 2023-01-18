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
export ParticleSoA # Particles represented as struct of arrays
export specieTable # Maping specie to mass and charge
export reset!      # Resets particle positions to zero (except initial position)
export setinitpos! # Sets the initial position of particles
export setinitvel! # Sets the initial velocity of particles
export kineticenergy # Computes the non-rel. kinetic energy at all time steps

#-------------#   
# Structs     # 
#-------------#-----------------------------------------------------------------
mutable struct ParticleSoA
    pos    ::Array{wpFloat, 3}
    vel    ::Array{wpFloat, 3}
    species::Vector{wpInt}   # Particle specie identifier (e.g. electron, proton)
    alive  ::Vector{Bool}
    
    
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
        return new(positions, velocities, species, alive)
    end # constructor 
end # mutable struct ParticleSoA


#----------------------#
# Particle set-methods #
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
    _, npart, N = size(particles.pos) # Will fail for GCA-particles
    v = norm3(particles.vel) # Will fail for GCA-particles
    Ek = zeros(npart, N)
    for i = 1:npart
        mass = specieTable[particles.species[i], 1]
        @. Ek[i,:] = 0.5*mass*v[i,:]^2
    end # loop i
    return Ek
end # function kineticenergy

#------------------#
# Global variables #
#------------------#------------------------------------------------------------
#              mass charge
specieTable = [m_e  -e     # Electron
               m_p   e     # Proton
               1.0 3.0     # Unit mass and charge = 3
               1.0 1.0]    # Unit mass and charge
end # module particles
