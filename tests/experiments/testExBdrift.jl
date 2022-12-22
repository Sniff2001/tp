#-------------------------------------------------------------------------------
# Created 19.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 particleDrift.jl
#
#-------------------------------------------------------------------------------
# Experiment testing if the ExB-drift of particles in a static homogeneous 
# electromagnetic field follows the analytical solution.
#-------------------------------------------------------------------------------

# Import external libraries
#import Pkg; Pkg.add("Plots")
#using Plots
using LinearAlgebra
using Test

# Import internal modules from tp/src
using Meshes
using Patches
using Particles
using Solvers
using Schemes
using Interpolations
using WorkingPrecision: wpFloat, wpInt

#-------------------------------------------------------------------------------
# EXPERIMENT PARAMTERS
#
numDims = 3       # Number of spatial dimensions
numParticles = 2  # Number of particles to simulae
species = [1, 2]  # Specifies the species of the particles 
                  #   (electron, proton)
tf = 2*2π/(5.93096958e7) # s. End time of simulation
                  #   multiple of electron (1eV) gyrofreq. at Bz=3.3e-4 G
numSteps = 1000   # Number of timesteps in the simulation
dt = tf/numSteps  # Size of timestep

#-------------------------------------------------------------------------------
# MESH CREATION
# The electromagnetic field is static and homogeneous. Hence the mesh only needs
# one cell, i.e. the boundary.
Bx = 0.0  # magnetic field component only in the positive z-direction
By = 0.0
Bz = 3.37213e-4 # G (so that electron larmor radius is 1 cm)
Ex = 0.0
Ey = 0.02Bz/tf  # Electric field component only in the positive y-direction
Ez = 0.0
B = zeros(wpFloat, numDims, 2, 2, 2)
E = zeros(wpFloat, numDims, 2, 2, 2)
B[3, :, :, :] .= Bz
E[1, :, :, :] .= Ex
E[2, :, :, :] .= Ey
E[3, :, :, :] .= Ez
# Create coordinates of cell
xCoords = [0.0, 1.0]
yCoords = [0.0, 1.0]
zCoords = [0.0, 1.0]
# Create Mesh instance
mesh = Mesh(B, E, xCoords, yCoords, zCoords)

#-------------------------------------------------------------------------------
# PARTICLE CONDITIONS
# Set initial position and velocities
vdriftx = Ey/Bz
# Electron initial velocity
vxe = 0.0#593_096.95848 # m/s (1eV electron
vye = 0.0
vze = 0.0
# Proton initial velocity
vxp = 13_841.122177 # m/s (1eV proton)
vyp = 0.0
vzp = 0.0
# Electron initial position
x0e = 0.0
y0e = 0.0
z0e = 0.0
# Proton initial position
x0p = 0.0
y0p = 0.0
z0p = 0.0
# Create array to contain initial conditions
pos = zeros(Float64, numDims, numParticles) # Position -> origin
vel = zeros(Float64, numDims, numParticles) # Velocity
vel[:, 1] = [vxe, vye, vze]  # Initial velocity electron
vel[:, 2] = [vxp, vyp, vzp]  # Initial velocity proton
pos[:, 1] = [x0e, y0e, z0e]  # Initial position electron
pos[:, 2] = [x0p, y0p, z0p]  # Initial position proton
# Create ParticlesSoA-instance
particles = ParticleSoA(pos, vel, species, numSteps)

#-------------------------------------------------------------------------------
# CREATE PATCH
patch = Patch(mesh,
              particles,
              Solvers.fullOrbit,
              Schemes.eulerCromer,
              Interpolations.trilinear,
              dt,
              numSteps,
              numParticles)

#-------------------------------------------------------------------------------
# RUN SIMULATION
Patches.run!(patch)

#-------------------------------------------------------------------------------
# ANALYTICAL SOLUTION 
times = LinRange(0, tf, numSteps + 1)

# Gyrofrequency
function ω(q, B, m)
    return q*B/m
end
# Larmor radius
function larmorRadius(ω, vperp)
    return vperp/ω
end
# x-position
function x(t, x0, ωt, rL, Ey, Bz)
    return @. rL*sin(ωt*t) + Ey/Bz*t
end
# y-positoin
function y(t, y0, ωt, rL, Ex, Bz)
    return @. -rL + rL*cos(ωt*t) - Ex/Bz*t
end
# z-position
function z(t, z0, vz0, ω, Ez, Bz)
    return @. z0 + 0.5ω/Bz*Ez*t^2 + vz0*t
end

# Define gyrofrequency and Larmor radius for both particles
qe = specieTable[1,2] # Electron charge
qp = specieTable[2,2] # Proton charge
me = specieTable[1,1] # Electron mass
mp = specieTable[2,1] # Proton mass
ωe = ω(qe, Bz, me)
ωp = ω(qp, Bz, mp)
vperpe = vxe - vdriftx # Subtract drift velocity to find vperp
vperpp = vxp - vdriftx
rLe = larmorRadius(ωe, vperpe)
rLp = larmorRadius(ωp, vperpp)

# Compute time evolution of position of both particles
pose = zeros(Float64, numDims, numSteps + 1)
posp = zeros(Float64, numDims, numSteps + 1)
pose[1, :] .= x(times, x0e, ωe, rLe, Ey, Bz)
pose[2, :] .= y(times, y0e, ωe, rLe, Ex, Bz)
pose[3, :] .= z(times, z0e, vze,  ωe, Ez, Bz)
posp[1, :] .= x(times, x0p, ωp, rLp, Ey, Bz)
posp[2, :] .= y(times, y0p, ωp, rLp, Ex, Bz)
posp[3, :] .= z(times, z0p, vzp, ωp, Ez, Bz)

#-------------------------------------------------------------------------------
# CALCULATE DEVIATIONS FROM ANALYTICAL SOLUTION
# Calculate the root mean squared error between numerical and analytical 
#   trajectories for both electron and proton for all position components.
# Electron
numPose = patch.tp.pos[:, 1, :]
numPosp = patch.tp.pos[:, 2, :]
rmsErrex = √(sum((pose[1,2:numSteps] .- numPose[1,2:numSteps]).^2)/numSteps)
rmsErrey = √(sum((pose[2,2:numSteps] .- numPose[2,2:numSteps]).^2)/numSteps)
rmsErrez = √(sum((pose[3,2:numSteps] .- numPose[3,2:numSteps]).^2)/numSteps)
# Proton
rmsErrpx = √(sum((posp[1,2:numSteps] .- numPosp[1,2:numSteps]).^2)/numSteps)
rmsErrpy = √(sum((posp[2,2:numSteps] .- numPosp[2,2:numSteps]).^2)/numSteps)
rmsErrpz = √(sum((posp[3,2:numSteps] .- numPosp[3,2:numSteps]).^2)/numSteps)

#-------------------------------------------------------------------------------
# TEST RESULTS
@testset verbose = true "ExB-drift" begin
    @test rmsErrex == 6.066588298648184e-5
    @test rmsErrey == 5.590378934365529e-5
    @test rmsErrez == 0.0
    @test rmsErrpx == 7.641789560048872e-14
    @test rmsErrpy == 3.369407609089028e-8
    @test rmsErrpz == 0.0
end # testset ExB-drift

#-------------------------------------------------------------------------------
# PLOT RESULTS
#p1 = plot(patch.tp.pos[1, 1, :], patch.tp.pos[2, 1, :], label="Numerical",
#          title="Electron", xlabel="x, m", ylabel="y, m")
#p1 = plot!(pose[1, :], pose[2, :], label="Analytical", ls=:dash)
#p2 = plot(patch.tp.pos[1, 2, :], patch.tp.pos[2, 2, :], label="Numerical",
#          title="Proton", xlabel="x, m", ylabel="y, m")
#p2 = plot!(posp[1, :], posp[2, :], label="Analytical", ls=:dash)
#
#plot(p1, p2, layout=(2,1))
