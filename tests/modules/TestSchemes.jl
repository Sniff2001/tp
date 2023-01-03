#-------------------------------------------------------------------------------
# Created 02.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#             TestSchemes.jl
#
#-------------------------------------------------------------------------------
# Module containing test of the Schemes-module.
#-------------------------------------------------------------------------------

module TestSchemes

using Test

using WorkingPrecision: wpInt, wpFloat
using Schemes

"""
"""

export testeuler
export testeulerCromer

# Test parameters
pos = [1, 2, 2, 2.5]
vel = [1, 2, 3, 1.5]
acc = [1, 3, 4, 10.5]
dt = 0.1

function testeuler(verbose::Bool)
    @testset verbose=verbose "euler" begin
        posAnswer = [1.1, 2.2, 2.3, 2.65]
        velAnswer = [1.1, 2.3, 3.4, 2.55]
        posNext, velNext = Schemes.euler(pos, vel, acc, dt)
        @test posAnswer == posNext
        @test velAnswer == velNext
    end # testset euler
end # function testeuler


function testeulerCromer(verbose::Bool)
    @testset verbose=verbose "eulerCromer" begin
        posAnswer = [1.11, 2.23, 2.34, 2.755]
        velAnswer = [1.1,   2.3,  3.4,  2.55]
        posNext, velNext = Schemes.eulerCromer(pos, vel, acc, dt)
        @test posAnswer == posNext
        @test velAnswer == velNext
    end # testset euler
end # function testeuler

end # module TestSchemes
        

