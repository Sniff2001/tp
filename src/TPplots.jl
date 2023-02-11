#-------------------------------------------------------------------------------
# Created 19.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                TPplots.jl
#
#-------------------------------------------------------------------------------
# Module containing plotting functions
#-------------------------------------------------------------------------------

module TPplots

# External libraries
import Plots
import PyPlot
# Internal libraries
using WorkingPrecision: wpInt, wpFloat
using Meshes
using Particles
using Utilities

export plotenergydistr
export plotplasmoid

#------#
# Mesh #
#-------------------------------------------------------------------------------
function quiverslice(mesh::Mesh, 
                     normal::String, 
                     point::wpInt,
                     field::String="B")
    if field == "B"
        f = mesh.bField
    elseif field == "E"
        f = mesh.eField
    else
        println("Error: Invalid field.")
    end
    if normal == "x"
        PyPlot.quiver(mesh.yCoords,
               mesh.zCoords,
               transpose(f[2, point, :, :]),
               transpose(f[3, point, :, :])
               )
        PyPlot.xlabel("y")
        PyPlot.ylabel("z")
    elseif normal == "y"
        PyPlot.quiver(mesh.xCoords,
               mesh.zCoords,
               transpose(f[1, :, point, :]),
               transpose(f[3, :, point, :])  
               )
        PyPlot.xlabel("x")
        PyPlot.ylabel("z")
    elseif normal == "z"
        PyPlot.quiver(mesh.xCoords,
               mesh.yCoords,
               transpose(f[1, :, :, point]),
               transpose(f[2, :, :, point])  
               )
        PyPlot.xlabel("x")
        PyPlot.ylabel("y")
    else
        println("Error: Plane not valid.")
    end
end # function quiverslice


function streamplotslice(mesh::Mesh, 
                         normal::String, 
                         point::wpInt,
                         field::String="B")
    if field == "B"
        f = mesh.bField
    elseif field == "E"
        f = mesh.eField
    else
        println("Error: Invalid field.")
    end
    if normal == "x"
        PyPlot.streamplot(mesh.yCoords,
               mesh.zCoords,
               transpose(f[2, point, :, :]),
               transpose(f[3, point, :, :])
               )
        xlabel("y")
        ylabel("z")
    elseif normal == "y"
        PyPlot.streamplot(mesh.xCoords,
               mesh.zCoords,
               transpose(f[1, :, point, :]),
               transpose(f[3, :, point, :])  
               )
        xlabel("x")
        ylabel("z")
    elseif normal == "z"
        PyPlot.streamplot(mesh.xCoords,
               mesh.yCoords,
               transpose(f[1, :, :, point]),
               transpose(f[2, :, :, point])  
               )
        xlabel("x")
        ylabel("y")
    else
        println("Error: Plane not valid.")
    end
end # function streamplotslice


#-----------#
# Particles #
#-------------------------------------------------------------------------------
function plotenergydistr(particles::ParticleSoA, 
                         snap     ::wpInt,
                         numbins  ::wpInt,
                         title
                         )
    absvel = norm2(particles.vel[:, :, snap])
    binrange = range(minimum(absvel), maximum(absvel), length=numbins)
    Plots.histogram(absvel, bins=binrange)
    Plots.xlabel!("Absolute velocity, m/s")
    Plots.xlabel!("Number of particles")
    Plots.title!(title)
end # function plotenergydistr
#
function plotenergydistr(absvel ::Vector{wpFloat},
                         numbins::wpInt,
                         title
                         )
    binrange = range(minimum(absvel), maximum(absvel), length=numbins)
    Plots.histogram(absvel, bins=binrange)
    Plots.xlabel!("Absolute velocity, m/s")
    Plots.xlabel!("Number of particles")
    Plots.title!(title)
end # function plotenergydistr

function plotplasmoid(pos,
                      vel,
                      times,
                      Ek,
                      xx,
                      yy, 
                      zz,
                      B,
                      numparticles
                      )
    for i = 1:numparticles
        PyPlot.plot(times, Ek[i, :])
    end
    PyPlot.title("Kinetic energy")
    PyPlot.xlabel("Time, s")
    PyPlot.ylabel("Energy, J")
    
    PyPlot.figure()
    PyPlot.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]))
    for i = 1:numparticles
        PyPlot.plot(pos[1,i,1], pos[2,i,1], marker="o")
    end
    PyPlot.quiver(pos[1, :, 1], pos[2, :, 1],
           vel[1, :, 1], vel[2, :, 1],
           width=0.003)
    PyPlot.title("Initial positions")
    
    PyPlot.figure()
    PyPlot.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]))
    for i = 1:numparticles
        PyPlot.scatter(pos[1,i,:], pos[2,i,:],
                s=0.1, marker=".")
    end
    PyPlot.title("Particle trajectories")
end # function plotplasmoid

end # module TPplots
