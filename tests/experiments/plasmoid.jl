#-------------------------------------------------------------------------------
# Created 13.01.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 plasmoid.jl
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Import external libraries
#import Pkg; Pkg.add("Plots")
#import Pkg; Pkg.add("PythonPlot")
#ENV["PYTHON"]="/mn/alruba2/astro/local/mamba/envs/py310/python"
#import Pkg; Pkg.add("PyPlot")
#ENV["PYTHON"]="python"
#Pkg.build("PyCall")
using PyPlot
#using Plots; pythonplot()
using LinearAlgebra

# Import internal modules from tp/src
using WorkingPrecision: wpFloat, wpInt
using Meshes
using Patches
using Particles
using Solvers
using Schemes
using Interpolations
using Utilities

"""
"""

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# SPATIAL PARAMETERS (x, y, z)
#
numdims = 3
# Lower bounds of the three spatial axes
xi0 = (0., 0., 0.)
# Upper bound of the three spatial axes
xif = (1., 1., 1.)
# Grid resolution of the axes
n = (10, 10, 2)

#-------------------------------------------------------------------------------
# SET THE MAGNETIC FIELD
#
bamp = 1.0
bconst = [10., 0., 0]
# Expectation and stdvalues. For creating vector potential with z-component
# normally distributed in x and y
x0, y0, z0 = xi0
xf, yf, zf = xif
μx = (xf - x0)/2 
μy = (yf - y0)/2 
σx = (xf - x0)/8
σy = (yf - y0)/8
μ  = (μx, μy)
σ  = (σx, σy)
amplitude = 10.0 # Amplitude of the z-component
axes, gridsizes, A = Utilities.normal3Donlyz(xi0, xif, n, μ, σ, amplitude)
# Derive the magnetic field from the curl of the vector-potential
B = Schemes.curl(A, gridsizes, Schemes.derivateCentral)
# Scale the field
@. B = bamp*B + bconst

#-------------------------------------------------------------------------------
# SET THE ELECTRIC FIELD
#
# Set the strength of the electric field components. Will be static and
# homogeneous 
Ex = 0.0
Ey = 0.0
Ez = 0.0

#-------------------------------------------------------------------------------
# PARTICLE CONDITIONS
#
numparticles = 20  # Number of particles to simulate
species = 4*ones(wpInt, numparticles)  # Specifies the species of the particles 
# Set initial position and velocities
# Initial velocity
vel0 = [0.5, 0.0, 0.0]
velf = [0.5, 0.0, 0.0]
# Initial position span
pos0 = [x0, y0, z0]
posf = [x0, yf, z0]
# Experiment conditions 
dt = 0.0001
tf = 2.0
numSteps = trunc(wpInt, tf/dt)   # Number of timesteps in the simulation
times = collect(0:dt:tf)
l = length(times)
println("Number of time steps = $numSteps. Length of times = $l")

#-------------------------------------------------------------------------------
# MESH CREATION
E = zeros(wpFloat, size(B))
E[1, :, :, :] .= Ex
E[2, :, :, :] .= Ey
E[3, :, :, :] .= Ez
# Create Mesh instance
xx, yy, zz = axes
mesh = Mesh(B, E, xx, yy, zz)

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
#pos = zeros(wpFloat, numdims, numParticles) # Position -> origin
#vel = zeros(wpFloat, numdims, numParticles) # Velocity
#vel[:, 1] = [vx0, vy0, vz0]  # Initial velocity electron
#pos[:, 1] = [posx0, posy0, posz0]  # Initial position electron

pos, vel = Utilities.initparticlesuniform(numparticles, pos0, posf, vel0, velf)
# Create ParticlesSoA-instance
particles = ParticleSoA(pos, vel, species, numSteps)

#-------------------------------------------------------------------------------
# CREATE PATCH
# Non-relativisitc Euler-Cromer
patch = Patch(mesh,
              particles,
              Solvers.fullOrbit,
              Schemes.eulerCromer,
              Interpolations.trilinear,
              dt,
              numSteps,
              numparticles)

#-------------------------------------------------------------------------------
# RUN SIMULATION
Patches.run!(patch)

#-------------------------------------------------------------------------------
# CALCULATE ENERGY
Ek = kineticenergy(particles)

#-------------------------------------------------------------------------------
# PLOT RESULTS
for i = 1:numparticles
    plot(times, Ek[i, :])
end
title("Kinetic energy")
xlabel("Time, s")
ylabel("Energy, J")

figure()
streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]))
for i = 1:numparticles
    plot(patch.tp.pos[1,i,1], patch.tp.pos[2,i,1], marker="o")
end
quiver(patch.tp.pos[1, :, 1], patch.tp.pos[2, :, 1],
       patch.tp.vel[1, :, 1], patch.tp.vel[2, :, 1],
       width=0.003)
title("Initial positions")

figure()
streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]))
for i = 1:numparticles
    scatter(patch.tp.pos[1,i,:], patch.tp.pos[2,i,:],
            s=0.1, marker=".")
end
title("Particle trajectories")



