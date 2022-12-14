#-------------------------------------------------------------------------------
#
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#
# Last edited: 12.12.22
#
#-------------------------------------------------------------------------------
#
#                 Patches.jl
#
#-------------------------------------------------------------------------------
# Module containing the Patch struct and methods
#-------------------------------------------------------------------------------

module Patches

using WorkingPrecision: wpFloat, wpInt
using Meshes
using Particles
using Solvers
using Interpolations

#-------------#   
# Export      # 
#-------------#-----------------------------------------------------------------
export Patch
export run

#-------------#   
# Main struct # 
#-------------#-----------------------------------------------------------------
mutable struct Patch
    mesh        ::Mesh
    tp          ::ParticleSoA # The trace particles
    solver      ::Function
    scheme      ::Function
    interpolator::Function
    dt          ::wpFloat
    numSteps    ::wpInt
    numParticles::wpInt
end

#---------#
# Methods #
#---------#---------------------------------------------------------------------
function run(patch::Patch)
    for i = 1:patch.numSteps # Over timesteps
        for j = 1:patch.numParticles # Over particles
            pos = patch.tp.pos[:, j]
            vel = patch.tp.vel[:, j]
            bField, eField = Interpolations.grid(patch.mesh,
                                                 patch.interpolator,
                                                 pos)
            acc = patch.solver(pos,
                               vel,
                               patch.tp.specie[j],
                               bField,
                               eField,
                               )
            pos, vel = patch.scheme(pos, vel, acc, patch.dt)
            patch.tp.pos[:, j] = pos
            patch.tp.vel[:, j] = vel
       end # loop over particles (j)
    end # loop over timesteps (i)
end # function: run
#-------------------------------------------------------------------------------

end # module
