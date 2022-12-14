#-------------------------------------------------------------------------------
#
# Created 13.12.22
# Author: e.s.oyre@astro.uio.no
#
# Last edited: 14.12.22
#
#-------------------------------------------------------------------------------
#
#             TestInterpolations.jl
#
#-------------------------------------------------------------------------------
# Module containing test of the Interpolations-module.
#-------------------------------------------------------------------------------

module TestInterpolations

using Test

using Interpolations
using Meshes

#export testlocateCell
export testtrilinear

function testtrilinear(verbose::Bool)
    @testset verbose=verbose "locateCell" begin
        #
        # Test interpolation of uniform field
        #
        B0 = ones(3,3,3,3) # uniform magnetic field
        E0 = ones(3,3,3,3) # uniform electric field
        mesh0 = Mesh(B0, E0) # Create mesh
        cellIndices00 = (1, 1, 1)
        cellIndices01 = (2, 2, 2)
        point00 = (0.0, 0.0, 0.0)
        point01 = (2/3, 2/3, 2/3)
        answer0 = ([1.0, 1.0, 1.0], [1.0, 1.0, 1.0])

        result00 = trilinear(mesh0, cellIndices00, point00)
        result01 = trilinear(mesh0, cellIndices01, point01)

        @test result00[1] ≈ answer0[1]
        @test result00[2] ≈ answer0[2]
        @test result01[1] ≈ answer0[1]
        @test result01[2] ≈ answer0[2]

        #
        # Test interpolation of linear field gradient
        #
        N = 3 # Number of vertices (grid points) on each axis of the mesh
        dx = 1/(N-1) # grid spacing

        # Field function
        linearField(x,y,z,a,b,c,d) = a*x + b*y + c*z .+ d
        a = b = c = d = 1 # Constants
        # Generate cartesian axes
        xx = collect(LinRange(0,1,N))
        yy = collect(LinRange(0,1,N))
        zz = collect(LinRange(0,1,N))
        # Evaluate fields on grid points
        B1 = zeros(3, N, N, N)
        E1 = zeros(3, N, N, N)
        # Broadcast into meshgrid?
        for i = 1:N
            for j = 1:N
                for k = 1:N
                    for dim = 1:3
                        B1[dim, i,j,k] = linearField(xx[i],
                                                     yy[j],
                                                     zz[k],
                                                     a, b, c, d)
                        E1[dim, i,j,k] = linearField(xx[i],
                                                     yy[j],
                                                     zz[k],
                                                     a, b, c, d)
                    end # loop over dim
                end # loop k
            end # loop j
        end # loop i

        # Initialise mesh
        mesh1 = Mesh(B1, E1, xx, yy, zz)
        # Define points to interpolate.
        #   Here the middle point of all cells.
        pointx = collect(0:N-1) .- dx/2
        pointy = collect(0:N-1) .- dx/2
        pointz = collect(0:N-1) .- dx/2
        # Evaluate the correct field values for these points
        answer = zeros(3, N-1, N-1, N-1) 
        # Loop through all cells of the mesh.
        #    Each index triple represent the corner of a celle
        for i = 1:N-1
            for j = 1:N-1
                for k = 1:N-1
                    answer = linearField(pointx[i],
                                         pointy[j],
                                         pointz[k],
                                         a, b, c, d)
                    # Interpolate
                    result = trilinear(mesh1, (i,j,k), (pointx[i],
                                                        pointy[j],
                                                        pointz[k]))
                    # Check result
                    @test result[1] ≈ [answer, answer, answer]
                    @test result[2] ≈ [answer, answer, answer]
                end # over k
            end # over j
        end # over i
    end # testset
end # function testtrilinear
end # module TestInterpolations
