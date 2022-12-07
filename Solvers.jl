#-------------------------------------------------------------------------------
#
# Created 02.12.22
# By Eilif @RoCS UiO.
# email: e.s.oyre@astro.uio.no
#
# Last edited: 07.12.22
#
#-------------------------------------------------------------------------------
#
#                 Solvers.jl
#
#-------------------------------------------------------------------------------
# Module containing the various solvers.
#-------------------------------------------------------------------------------
module Solvers

# Standard library
using LinearAlgebra: cross
# Internal modules
using Particles: specieTable


function fullOrbit(pos, vel, specie, bField, eField)
    mass   = specieTable[specie, 1]
    charge = specieTable[specie, 2]
    acc = charge/mass * (eField + cross(vel, bField))
end

end
