# Created 31.03.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 gradBdrift.jl
#
#-------------------------------------------------------------------------------
# This experiment tests both the full orbit and GCA solvers against gradient B
# drift.
#-------------------------------------------------------------------------------

using LinearAlgebra
using Test

# Import internal modules from tp/src
using WorkingPrecision: wpFloat, wpInt
using Constants:        k_B
using Meshes:           Mesh
using Patches:          Patch, run!
using Particles:        ParticleSoA, GCAParticleSoA, kineticenergy, specieTable
using Solvers          
using Schemes
using Interpolations
using Utilities

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
dt = 0.02          # Time step [s]                      |
tf = 5.0 #n=100   # End time of simulation [s]         |
#tf = 2.3 #n=10   # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
species = 4*ones(wpInt, numparticles) # Specifies the species of the particles 

#...............................................
# INITIAL CONDITIONS
pos0 = [0.5, 0.5, 0.5]
vel0 = [0.0, 0.3, 0.0]

#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
xi0 = (0., 0., 0.)
# Upper bound of the three spatial axes
xif = (1., 1., 1.)
# Grid resolution of the axes
#n = (10, 10, 2)
ni = (100, 100, 2)
 
#...............................................
# MAGNETIC FIELD PARAMETERS
a = 10.0 # gradient in magnetic field
B0 = 0.0 # Additional constant
function gradBfield(
    x::wpFloat,
    y::wpFloat,
    z::wpFloat
    )
    return [0.0, 0.0, a*y + B0]
end

#...............................................
# ELECTRIC FIELD PARAMETERS
# It will not exist
Ex = 0.0
Ey = 0.0

#...............................................
# SOLVER CONDITIONS
FOscheme = Schemes.rk4  # Scheme for integrating the full-orbit eqs.
GCAscheme = Schemes.euler # Scheme for integrating the GCA eqs.
FOinterp = Interpolations.trilinear # Interpolation scheme for full orbit
GCAinterp = Interpolations.trilinearGCA # Interpolation scheme for GCA
pbc    = (true, true, true) # (x,y,z) Are mesh boundary conditions periodic?


#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE AXES, MAGNETIC FIELD AND ELECTRIC FIELD
xx, yy, zz, dx, dy, dz = Utilities.createaxes(xi0, xif, ni)
Bfield = zeros(wpFloat, numdims, ni[1], ni[2], ni[3])
Efield = zeros(size(Bfield))
Efield[1,:,:,:] .= Ex
Efield[2,:,:,:] .= Ey
#Efield[1,:,1:30,:] .= 0.0
Utilities.discretise!(Bfield, xx, yy, zz, gradBfield)

#-------------------------------------------------------------------------------
# MESH CREATION
# Create Mesh instance
mesh = Mesh(Bfield, Efield, xx, yy, zz)

#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numsteps = trunc(wpInt, tf/dt)   # Number of timesteps in the simulation
println("Number of time steps = $numsteps.")

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
(B⃗, E⃗), _ = gridinterp(mesh,
                       trilinear,
                       pos0)
B = norm(B⃗)
b̂ = B⃗/B
v = norm(vel0)
vparal = [vel0 ⋅ b̂]
vperp = √(v^2 - vparal[1]^2)
mass = specieTable[species[1], 1]
μ = [mass*vperp^2/(2B)]

tpFO = ParticleSoA(pos0, vel0, species, numsteps)
tpGCA = GCAParticleSoA(pos0, vparal, μ, species, numsteps)

#-------------------------------------------------------------------------------
# CREATE PATCH
# Non-relativisitc Euler-Cromer
patchFO = Patch(mesh,
                tpFO,
                Solvers.fullOrbit,
                FOscheme,
                FOinterp,
                dt,
                numsteps,
                numparticles,
                pbc
                )
patchGCA = Patch(mesh,
                 tpGCA,
                 Solvers.GCA,
                 GCAscheme,
                 GCAinterp,
                 dt,
                 numsteps,
                 numparticles,
                 pbc
                 )

#...............................................................................
# RUN SIMULATION
@time run!(patchFO)
@time run!(patchGCA)
#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
charge = specieTable[species[1], 2]
vdrift = [-a*mass*vel0[2]^2/(2charge*B^2), 0.0, 0.0]
posf_anal = tf*vdrift + pos0
posf_num  = patchGCA.tp.R[:,1,end]
@testset verbose = true "GCA: Euler" begin
    @test posf_anal ≈ posf_num
end # testset GCA: Euler
    
