#-------------------------------------------------------------------------------
# Created 02.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#             TestPatches.jl
#
#-------------------------------------------------------------------------------
# Module containing test of the Patches-module.
#-------------------------------------------------------------------------------

module TestPatches

using Test

using WorkingPrecision: wpInt, wpFloat
using Patches

"""
There are currently no Patches-specific unit-tests.

See testExBdrift.jl for indirect testing of run!(patch::Patch).
"""
