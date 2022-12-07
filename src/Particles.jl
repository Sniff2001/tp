#-------------------------------------------------------------------------------
#
# Created 05.12.22
# By Eilif @RoCS UiO.
# email: e.s.oyre@astro.uio.no
#
# Last edited: 07.12.22
#
#-------------------------------------------------------------------------------
#
#                Particles.jl
#
#-------------------------------------------------------------------------------
# Module containing the Particle structs and methods
#-------------------------------------------------------------------------------

module Particles

using WorkingPrecision: wpFloat, wpInt


#-------------#   
# Exports     # 
#-------------#-----------------------------------------------------------------
export ParticleSoA # Particles represented as struct of arrays
export specieTable # Maping specie to mass and charge

#-------------#   
# Structs     # 
#-------------#-----------------------------------------------------------------
mutable struct ParticleSoA
    pos   ::Matrix{wpFloat} # Particle positions with shape (numDims x numPart)
    vel   ::Matrix{wpFloat} # Particles velocities    shape (numDims x numPart)
    specie::Vector{wpInt}   # Particle specie identifier (e.g. electron, proton)
end
#-------------------------------------------------------------------------------

#------------------#
# Global variables #
#------------------#------------------------------------------------------------
specieTable = [1.0  -1.0   # Electron
               10.0  1.0]  # Proton
end
