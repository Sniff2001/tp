#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
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
using Schemes
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
    periodicBC  ::Tuple{Bool, Bool, Bool}

    # Constructors
    #--------------------------------------------------------------------------
    function Patch(mesh        ::Mesh,
                   tp          ::ParticleSoA, # The trace particles
                   solver      ::Function,
                   scheme      ::Function,
                   interpolator::Function,
                   dt          ::wpFloat,
                   numSteps    ::wpInt,
                   numParticles::wpInt
                   )
        new(mesh, 
            tp, 
            solver, 
            scheme, 
            interpolator, 
            dt, 
            numSteps, 
            numParticles,
            (false, false, false)
            )
    end # constructor

    function Patch(mesh        ::Mesh,
                   tp          ::ParticleSoA, # The trace particles
                   solver      ::Function,
                   scheme      ::Function,
                   interpolator::Function,
                   dt          ::wpFloat,
                   numSteps    ::wpInt,
                   numParticles::wpInt,
                   periodicBC  ::Tuple{Bool, Bool, Bool}
                   )
        new(mesh, 
            tp, 
            solver, 
            scheme, 
            interpolator, 
            dt, 
            numSteps, 
            numParticles,
            periodicBC
            )
    end # constructor
end # mutable struct Patch

#---------#
# Methods #
#---------#---------------------------------------------------------------------
function run!(patch::Patch)
    for i = 1:patch.numSteps # Over timesteps
        for j = 1:patch.numParticles # Over particles
            if patch.tp.alive[j] == false
                continue
            end # if !alive
            pos = patch.tp.pos[:, j, i]
            vel = patch.tp.vel[:, j, i]
            pos, vel = patch.solver(pos,
                                    vel,
                                    patch.tp.species[j],
                                    patch.mesh,
                                    patch.dt,
                                    patch.interpolator,
                                    patch.scheme,
                                    )
            for k = 1:patch.mesh.numdims
                if pos[k] < patch.mesh.domain[k, 1]
                    if patch.periodicBC[k] == true
                        pos[k] = patch.mesh.domain[k, 2] +
                            (pos[k] - patch.mesh.domain[k, 1])
                    else
                        patch.tp.alive[j] = false # kill particle
                        patch.tp.pos[:, j, i+2:patch.numSteps+1] .= 
                            patch.tp.pos[:, j, i]
                        patch.tp.vel[:, j, i+2:patch.numSteps+1] .= 
                            patch.tp.vel[:, j, i]
                        pos = patch.tp.pos[:, j, i]
                        vel = patch.tp.vel[:, j, i]
                        break
                   end
                elseif pos[k] > patch.mesh.domain[k, 2]
                    if patch.periodicBC[k] == true
                        pos[k] = patch.mesh.domain[k, 1] +
                            (pos[k] - patch.mesh.domain[k, 2])
                    else
                        patch.tp.alive[j] = false # kill particle
                        patch.tp.pos[:, j, i+2:patch.numSteps+1] .= 
                            patch.tp.pos[:, j, i]
                        patch.tp.vel[:, j, i+2:patch.numSteps+1] .= 
                            patch.tp.vel[:, j, i]
                        pos = patch.tp.pos[:, j, i]
                        vel = patch.tp.vel[:, j, i]
                        break
                    end # if particle alive
                end # if particle outside domain
            end # loop dimensions
            patch.tp.pos[:, j, i+1] = pos
            patch.tp.vel[:, j, i+1] = vel
       end # loop over particles (j)
    end # loop over timesteps (i)
end # function: run
#-------------------------------------------------------------------------------

end # module
