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
using Interpolations_tp

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
    tp          ::TraceParticle # The trace particles
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
                   tp          ::TraceParticle, # The trace particles
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
                   tp          ::TraceParticle, # The trace particles
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
        Particles.push!(patch.tp,
                        patch.mesh,
                        i,
                        patch.dt,
                        patch.solver,
                        patch.interpolator,
                        patch.scheme,
                        patch.periodicBC,
                        )
    end # loop over timesteps (i)
end # function: run


#----------------#
# Base functions #
#-------------------------------------------------------------------------------
"""
    Base.copy(p::Patch)
Make a deep copy of a Patch-type.
"""
function Base.copy(p::Patch)
    Patch(p.mesh,
          copy(p.tp),
          p.solver,
          p.scheme,
          p.interpolator,
          p.dt,
          p.numSteps,
          p.numParticles,
          p.periodicBC
          )
end # function Base.copy


function Base.Multimedia.display(p::Patch)
    println("""Instance of mutable struct:
    Patches.Patch
        mesh        ::Meshes.Mesh
        tp          ::Particles.TraceParticle 
        solver      ::Function
        scheme      ::Function
        interpolator::Function
        dt          ::wpFloat
        numSteps    ::wpInt
        numParticles::wpInt""")
end # function Base.Multimedia.display
#-------------------------------------------------------------------------------

end # module Patches
