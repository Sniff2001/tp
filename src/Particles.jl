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

using LinearAlgebra:    norm, ⋅

using WorkingPrecision: wpFloat, wpInt
using Constants:        m_e, m_p, e
using Utilities:        norm3
using Meshes
using Interpolations_tp


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
export computeμ

#------------------#
# Global variables #
#------------------#------------------------------------------------------------
#              mass charge
specieTable = [m_e   -e     # Electron
               m_p    e     # Proton
               wpFloat(1.0)  wpFloat( 3.0)     # Unit mass and charge = 3
               wpFloat(1.0)  wpFloat( 1.0)    # Unit mass and charge
               wpFloat(1.0)  wpFloat(-1.0)    # Unit mass and negative unit charge
               ]

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
    function ParticleSoA( # "Default" constructor
        pos    ::Array{wpFloat, 3},
        vel    ::Array{wpFloat, 3},
        species::Vector{wpInt},   # Particle specie ID (e.g. electron, proton)
        alive  ::Vector{Bool},
        weight ::Vector{wpFloat}
        )
        return new(pos, vel, species, alive, weight)
    end # constructor 
    #|
    function ParticleSoA(
        pos0     ::Matrix{wpFloat},
        vel0     ::Matrix{wpFloat},
        species ::Vector{wpInt},
        numSteps::Integer
        )
        numDims, numParticles = size(pos0)
        positions  = zeros(wpFloat, numDims, numParticles, numSteps + 1)
        velocities = zeros(wpFloat, numDims, numParticles, numSteps + 1)
        positions[:, :, 1] .= pos0
        velocities[:, :, 1] .= vel0
        alive = ones(Bool, numParticles)
        weight = ones(wpFloat, numParticles)
        return new(positions, velocities, species, alive, weight)
    end # constructor 
    #|
    function ParticleSoA(
        pos0    ::Vector{wpFloat},
        vel0    ::Vector{wpFloat},
        species ::Vector{wpInt},
        numSteps::Integer
        )
        numParticles = 1
        numDims = length(pos0)
        positions  = zeros(wpFloat, numDims, numParticles, numSteps + 1)
        velocities = zeros(wpFloat, numDims, numParticles, numSteps + 1)
        positions[:, :, 1] .= pos0
        velocities[:, :, 1] .= vel0
        alive = ones(Bool, numParticles)
        weight = ones(wpFloat, numParticles)
        return new(positions, velocities, species, alive, weight)
    end # constructor 
end # mutable struct ParticleSoA


mutable struct GCAParticleSoA <: TraceParticle
    R      ::Array{wpFloat, 3} # The position of the guiding centre
    vparal ::Matrix{wpFloat}   # Velocity parallel to the magnetic field
    μ      ::Vector{wpFloat}   # Magnetic moment μ of particle
    species::Vector{wpInt}     # Particle specie ID (e.g. electron, proton)
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
    function GCAParticleSoA( # "Default" constructor
        R      ::Array{wpFloat, 3}, # The position of the guiding centre
        vparal ::Matrix{wpFloat},   # The velocity parallel to the magnetic field
        μ      ::Vector{wpFloat},   # Magnetic moment μ of particle
        species::Vector{wpInt},     # Particle specie identifier
                                    #   (e.g. electron, proton)
        alive  ::Vector{Bool},
        weight ::Vector{wpFloat}
        )
        return new(R, vparal, μ, species, alive, weight)
    end # constructor 
    #|
    function GCAParticleSoA( # Given only initial conditions
        R0      ::Matrix{wpFloat},
        vparal0 ::Vector{wpFloat},
        μ       ::Vector{wpFloat},
        species ::Vector{wpInt},
        numSteps::Integer
        )
        numDims, numParticles = size(R0)
        R      = zeros(wpFloat, numDims, numParticles, numSteps + 1)
        vparal = zeros(wpFloat, numParticles, numSteps + 1)
        R[:, :, 1] .= R0
        vparal[:, :, 1] .= vparal0
        alive = ones(Bool, numParticles)
        weight = ones(wpFloat, numParticles)
        return new(R, vparal, μ, species, alive, weight)
    end # constructor 
    #|
    function GCAParticleSoA( # Given only initial conditions
        R0      ::Vector{wpFloat},
        vparal0 ::Vector{wpFloat},
        μ       ::Vector{wpFloat},
        species ::Vector{wpInt},
        numSteps::Integer
        )
        numParticles = 1
        numDims = length(R0)
        R      = zeros(wpFloat, numDims, numParticles, numSteps + 1)
        vparal = zeros(wpFloat, numParticles, numSteps + 1)
        R[:, :, 1] .= R0
        vparal[:, :, 1] .= vparal0
        alive = ones(Bool, numParticles)
        weight = ones(wpFloat, numParticles)
        return new(R, vparal, μ, species, alive, weight)
    end # constructor
end # mutable struct ParticleSoA


#----------------#
# Base functions #
#-------------------------------------------------------------------------------
"""
    Base.copy(tp::TraceParticle)
Make a deep copy of a TraceParticle-type.
"""
function Base.copy(tp::ParticleSoA)
    ParticleSoA(copy(tp.pos),
                copy(tp.vel),
                copy(tp.species),
                copy(tp.alive),
                copy(tp.weight)
                )
end # function Base.copy
#|
function Base.copy(tp::GCAParticleSoA)
    GCAParticleSoA(copy(tp.R),
                   copy(tp.vparal),
                   copy(tp.μ),
                   copy(tp.species),
                   copy(tp.alive),
                   copy(tp.weight)
                   )
end # function Base.copy
                     

function Base.Multimedia.display(tp::ParticleSoA)
    println("""Instance of mutable struct:
    Particles.ParticleSoA <: Particles.TraceParticle
        pos    ::Array{wpFloat, 3}
        vel    ::Array{wpFloat, 3}
        species::Vector{wpInt}   
        alive  ::Vector{Bool}
        weight ::Vector{wpFloat}
    """)
    println("................................................................")
    println("GCAPaticleSoA.pos:")
    display(tp.pos)
    println("................................................................")
    println("GCAPaticleSoA.vel:")
    display(tp.vel)
    println("................................................................")
    println("GCAPaticleSoA.species:")
    display(tp.species)
    println("................................................................")
    println("GCAPaticleSoA.alive:")
    display(tp.alive)
    println("................................................................")
    println("GCAPaticleSoA.weight:")
    display(tp.weight)
end # function Base.Multimedia.display
#|
function Base.Multimedia.display(tp::GCAParticleSoA)
    println("""Instance of mutable struct:
    Particles.ParticleSoA <: Particles.TraceParticle
        R      ::Array{wpFloat, 3} 
        vparal ::Matrix{wpFloat}   
        μ      ::Vector{wpFloat}  
        species::Vector{wpInt}   
        alive  ::Vector{Bool}
        weight ::Vector{wpFloat}
    """)
    println("................................................................")
    println("GCAPaticleSoA.R:")
    display(tp.R)
    println("................................................................")
    println("GCAPaticleSoA.vparal:")
    display(tp.vparal)
    println("................................................................")
    println("GCAPaticleSoA.μ:")
    display(tp.μ)
    println("................................................................")
    println("GCAPaticleSoA.species:")
    display(tp.species)
    println("................................................................")
    println("GCAPaticleSoA.alive:")
    display(tp.alive)
    println("................................................................")
    println("GCAPaticleSoA.weight:")
    display(tp.weight)
end # function Base.Multimedia.display


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


#---------------------------#
# Push particles a timestep #
#-------------------------------------------------------------------------------
function push!(
    tp    ::GCAParticleSoA,
    mesh  ::Mesh,
    time  ::wpInt,
    dt    ::wpFloat,
    solver::Function,
    interp::Function,
    scheme::Function,
    periodicBC::Tuple{Bool, Bool, Bool},
    )
    for j in eachindex(tp.alive)
        if tp.alive[j] == false
            continue
        end
        pos = tp.R[:,j,time]
        vel = tp.vparal[j,time]
        pos, vel = solver(
            pos,
            vel,
            tp.μ[j],
            tp.species[j],
            mesh,
            dt,
            interp,
            scheme,
        )
        checkboundary!(pos, tp.alive[j], periodicBC, mesh.domain)
        tp.R[:,j,time+1] = pos
        tp.vparal[j,time+1] = vel
    end # loop over particles
end # function push!
#|
function push!(
    tp    ::ParticleSoA,
    mesh  ::Mesh,
    time  ::wpInt,
    dt    ::wpFloat,
    solver::Function,
    interp::Function,
    scheme::Function,
    periodicBC::Tuple{Bool, Bool, Bool},
    )
    for j in eachindex(tp.alive)
        if tp.alive[j] == false
            continue
        end
        pos = tp.pos[:,j,time]
        vel = tp.vel[:,j,time]
        pos, vel = solver(
            pos,
            vel,
            tp.species[j],
            mesh,
            dt,
            interp,
            scheme,
        )
        checkboundary!(pos, tp.alive[j], periodicBC, mesh.domain)
        tp.pos[:,j,time+1] = pos
        tp.vel[:,j,time+1] = vel
    end # loop over particles
end # function push!

function checkboundary!(
    pos       ::Vector{wpFloat},
    alive     ::Bool,
    periodicBC::Tuple{Bool, Bool, Bool},
    domain    ::Matrix{wpFloat},
    )
    for k in 1:length(domain[:,1])
        if pos[k] < domain[k,1]
            if periodicBC[k] == true
                pos[k] = domain[k, 2] + (pos[k] - domain[k, 1])
            else
                alive = false # kill particle
                break
            end
        elseif pos[k] > domain[k, 2]
            if periodicBC[k] == true
                pos[k] = domain[k, 1] + (pos[k] - domain[k, 2])
            else
                alive = false # kill particle
                break
            end # if particle alive
        end # if particle outside domain
    end
end # function checkboundary!

#----------------------#
# Auxiliary quantities #
#-------------------------------------------------------------------------------
function kineticenergy(particles::ParticleSoA)
    _, npart, N = size(particles.pos) 
    v = norm3(particles.vel) 
    Ek = zeros(npart, N)
    for j = 1:npart
        mass = specieTable[particles.species[j], 1]
        @. Ek[j,:] = 0.5*mass*v[j,:]^2
    end # loop j
    return Ek
end # function kineticenergy
#|
function kineticenergy(
    particles::ParticleSoA,
    mesh     ::Mesh,
    interp   ::Function
    )
    return kineticenergy(particles)
end 
#|
function kineticenergy(
    particles::GCAParticleSoA,
    mesh     ::Mesh,
    interp   ::Function
    )
    npart, N = size(particles.vparal) 
    mass = specieTable[particles.species[:], 1]
    vperp = getvperp(particles.R, particles.μ, mass, mesh, interp)
    Ek = zeros(npart, N)
    Ekparal = zeros(npart, N)
    Ekperp = zeros(npart, N)
    for j = 1:npart
        v2 = particles.vparal[j,:].^2 + vperp[j,:].^2
        @. Ek[j,:] = 0.5*mass[j]*v2
        @. Ekparal[j,:] = 0.5*mass[j]*particles.vparal[j,:]^2
        @. Ekperp[j,:] = 0.5*mass[j]*vperp[j,:]^2
    end # loop j
    return Ek, Ekparal, Ekperp
end # function kineticenergy


function getvperp(
    R     ::Array{wpFloat, 3},
    μ     ::Vector{wpFloat},
    mass  ::Vector{wpFloat},
    mesh  ::Mesh,
    interp::Function,
    )
    _, numparticles, numsteps = size(R)
    vperp = zeros(numparticles, numsteps)
    for i = 1:numsteps
        for j = 1:numparticles
            bfield, _ = gridinterp(mesh.bField, interp, R[:,j,i],
                                   mesh.xCoords, mesh.yCoords, mesh.zCoords
                                   )
            B = norm(bfield)
            vperp[j,i] = √(2μ[j]*B/mass[j])
        end # loop over j: particles
    end # loop over i: timesteps
    return vperp
end # function getvperp

"""
    computeμ(
        tp::ParticleSoA,
        mesh::Mesh,
        intep::Function,
        )
Calculates and returns μ for all particles at all time steps by
interpolating to the particle positions. This is a post-processing procudureand
is only accurate if the same interpolation method is used as in the simulation
itself.
"""
function computeμ(
    tp::ParticleSoA,
    mesh::Mesh,
    interp::Function,
    )
    m = specieTable[tp.species, 1]
    _, numparticles, numsteps = size(tp.pos)
    μ = zeros(numparticles, numsteps)
    for i = 1:numsteps
        for j = 1:numparticles
            v⃗ = tp.vel[:,j,i]
            r⃗ = tp.pos[:,j,i]
            bfield, _ = gridinterp(mesh.bField, interp, r⃗,
                                   mesh.xCoords, mesh.yCoords, mesh.zCoords
                                   )
            B = norm(bfield)
            b̂ = bfield/B
            vparal = v⃗⋅b̂
            vperp =  norm(v⃗ - vparal*b̂)
            μ[j,i] = 0.5*m[j]*vperp^2/B
        end # loop over j: particles
    end # loop over i: timesteps
    return μ
end # function computeμ

end # module particles
