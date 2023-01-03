#-------------------------------------------------------------------------------
# Created 02.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#             TestParticles.jl
#
#-------------------------------------------------------------------------------
# Module containing test of the Particles-module.
#-------------------------------------------------------------------------------

module TestParticles

using Test

using WorkingPrecision: wpInt, wpFloat
using Particles

"""
There are currently no Particles-specific unit-tests.

See testExBdrift.jl for indirect testing of 
    ParticleSoA(pos, vel species, numsteps)
"""
