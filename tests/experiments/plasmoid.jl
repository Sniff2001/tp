#-------------------------------------------------------------------------------
# Created 13.01.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 plasmoid.jl
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Import external libraries
#import Pkg; Pkg.add("Plots")
#import Pkg; Pkg.add("PythonPlot")
#import Pkg; Pkg.add("PyPlot")
using PyPlot
#using Plots; pythonplot()
using LinearAlgebra

# Import internal modules from tp/src
using WorkingPrecision: wpFloat, wpInt
using Meshes
using Patches
using Particles
using Solvers
using Schemes
using Interpolations

"""
"""

function normaldistr(x, μ=0, σ=√2)
    return @.  1/(σ*√2π)*exp(-0.5((x-μ)/σ)^2)
end # normaldistr


x0 = 0.
y0 = 0.
z0 = 0.

xf = 1.
yf = 1.
zf = 1.

nx = 10
ny = 10
nz = 2
ndims = 3

xx = collect(LinRange(x0, xf, nx))
dx = xx[2] - xx[1]
yy = collect(LinRange(y0, yf, ny))
dy = yy[2] - yy[1]
if nz == 1
    zz = [z0]
    dz = 0.
else
    zz = collect(LinRange(z0, zf, nz))
    dz = zz[2] - zz[1]
end

A = zeros(ndims, nx, ny, nz)
μx = (xf - x0)/2 
μy = (yf - y0)/2 
μz = (zf - z0)/2
σx = (xf - x0)/9
σy = (yf - y0)/9
σz = (zf - z0)/9

for i = 1:nx
    for j = 1:ny
        A[3,i,j,:] .= normaldistr(xx[i], μx, σx) * normaldistr(yy[j], μy, σy)
    end
end

#for i = 1:nx
#    for k = 1:nz
#        A[2,i,:,k] = normaldistr(yy, μy, σy)
#    end 
#end
#for j = 1:ny
#    for k = 1:nz
#        A[1,:,j,k] = normaldistr(xx, μx, σx)
#    end
#end

B = Schemes.curl(A, dx, dy, dz, Schemes.derivateCentral)
#@. B[1,:,:,:] = B[1,:,:,:] + 1
streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]))
