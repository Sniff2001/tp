
#-------------------------------------------------------------------------------
#
# Created 13.12.22
# Author: e.s.oyre@astro.uio.no
#
# Last edited: 13.12.22
#
#-------------------------------------------------------------------------------
# Script for running tests.
#-------------------------------------------------------------------------------

using Test
using TestInterpolations

@testset verbose = true "Interpolations" begin
    testtrilinear(true)
    testlocateCell(true)
end
