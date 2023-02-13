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
using LinearAlgebra:    norm
# Internal libraries
using WorkingPrecision: wpInt, wpFloat
using Meshes
using Patches:          Patch
using Particles
using Utilities:        rejectionsampling, norm4, norm2
using Interpolations
using Schemes
using Solvers

export plot

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
function plotenergydistr(particles::Particles.ParticleSoA, 
                         snap     ::wpInt,
                         numbins  ::wpInt,
                         title
                         )
    absvel = norm2(particles.vel[:, :, snap])
    binrange = range(minimum(absvel), maximum(absvel), length=numbins)
    Plots.histogram(absvel, bins=binrange)
    Plots.xlabel!("Absolute velocity, m/s")
    Plots.ylabel!("Number of particles")
    Plots.title!(title)
end # function plotenergydistr
#
function plotenergydistr(particles::Particles.GCAParticleSoA, 
                         snap     ::wpInt,
                         numbins  ::wpInt,
                         title
                         )
    vparal = particles.vparal[:, snap]
    binrange = range(minimum(vparal), maximum(vparal), length=numbins)
    Plots.histogram(vparal, bins=binrange)
    Plots.xlabel!("Parallel velcity, m/s")
    Plots.ylabel!("Number of particles")
    Plots.title!(title)
end # function plotenergydistr
function plotenergydistr(absvel ::Vector{wpFloat},
                         numbins::wpInt,
                         title
                         )
    binrange = range(minimum(absvel), maximum(absvel), length=numbins)
    Plots.histogram(absvel, bins=binrange)
    Plots.xlabel!("Absolute velocity, m/s")
    Plots.ylabel!("Number of particles")
    Plots.title!(title)
end # function plotenergydistr

function plotplasmoid(
    tp::TraceParticle,
    times,
    Ek,
    xx,
    yy, 
    zz,
    B,
    numparticles,
    labelling
    )

    pos = getpos(tp)
    # Energy plot
    PyPlot.figure()
    for i = 1:numparticles
        PyPlot.plot(times, Ek[i, :], label="#$i")
    end
    PyPlot.title("Kinetic energy")
    PyPlot.xlabel("Time, s")
    PyPlot.ylabel("Energy, J")
    if labelling
        PyPlot.legend()
    end
    
    # Streamplot of magnetic field with initial positions and velocity
    PyPlot.figure()
    PyPlot.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]))
    for i = 1:numparticles
        PyPlot.plot(pos[1,i,1], pos[2,i,1], marker="o")
    end
    if typeof(tp) == ParticleSoA
        # Plot starting velocity as an arrow
        vel = getvel(tp)
        PyPlot.quiver(pos[1, :, 1], pos[2, :, 1],
                      vel[1, :, 1], vel[2, :, 1],
                      width=0.003)
        PyPlot.title("Particle trajectories")
    end
    
    # Streamplot of magnetic field with particle trajectories
    PyPlot.figure()
    PyPlot.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]))

    cm = PyPlot.get_cmap(:tab20)
    colorrange = (0:numparticles) ./ numparticles

    for i = 1:numparticles
        plotperiodictrajectory(pos[:,i,:], i, cm, colorrange)
    end
    for i = 1:numparticles
        PyPlot.plot(pos[1,i,1], pos[2,i,1], marker="o", color="blue")
        # Mark position after t=1.5
        if i == 1
            PyPlot.plot(pos[1,1,1501], pos[2,1,1501], marker="^",
                        color="black",label="t = 1.5s", linestyle="None") 
        else
            PyPlot.plot(pos[1,i,1501], pos[2,i,1501], marker="^",
                        color="black") 
        end
    end

    if typeof(tp) == ParticleSoA
        # Plot starting velocity as an arrow
        vel = getvel(tp)
        PyPlot.quiver(pos[1, :, 1], pos[2, :, 1],
                      vel[1, :, 1], vel[2, :, 1],
                      width=0.003)
        PyPlot.title("Particle trajectories")
    end

    if labelling
        PyPlot.legend()
    end

end # function plotplasmoid


function plotperiodictrajectory(pos, partidx, cm, colorrange)
    extent = [1.0, 1.0, 1.0] # Hardcoded for now
    posjumps = diff(pos, dims=2)
    posjumps = @. abs(posjumps)
    mask = @. isapprox(posjumps, extent, rtol=0.01)
    indices = findall(mask[1:2, :])
    numindices = length(indices)
    if numindices > 0
        j = 1
        for s = 1:length(indices)
            i = indices[s][2] # We only care about what time step the particle
            # hit the boundary, not which boundary was crossed.
            PyPlot.plot(pos[1,j:i], pos[2,j:i], color=cm(colorrange[partidx]))
                        
            j = i+1
        end
        PyPlot.plot(pos[1,j:end], pos[2,j:end], color=cm(colorrange[partidx]),
                    label="#$partidx")
    else
        PyPlot.plot(pos[1,:], pos[2,:], color=cm(colorrange[partidx]),
                    label="#$partidx")
    end
end 


#---------#
# Patches #
#-------------------------------------------------------------------------------
function plot(
    patch    ::Patch,
    labelling::Bool=false
    )
    times = collect(range(0.0, step=patch.dt, length=patch.numSteps+1))
    Ek = kineticenergy(patch.tp)
    plotplasmoid(
        patch.tp,
        times,
        Ek,
        patch.mesh.xCoords, patch.mesh.yCoords, patch.mesh.zCoords,
        patch.mesh.bField,
        patch.numParticles,
        labelling
    )

    numbins = 20
    p1 = plotenergydistr(patch.tp, 1, numbins, "Initial energy")
    p2 = plotenergydistr(patch.tp, patch.numSteps+1, numbins, "final energy")
    #p3 = plotenergydistr(patch.tp.vel[1, :, 1], numbins, "Vx-initial")
    #p4 = plotenergydistr(patch.tp.vel[1, :, end], numbins, "Vx-final")
    Plots.plot(p1,p2)#p3,p4)
end


#---------#
# Other   #
#-------------------------------------------------------------------------------
function generatefieldline(
    mesh,
    initpos,
    stepsize,
    numsteps,
    interpolator,
    scheme
    )

    # Find out how many steps are needed in both directions
#    pos = zeros(mesh.numdims)
#    # Forward
#    pos .= initpos
#    stepsforward = -1
#    # While the position is till in the domain
#    while all(pos .>= mesh.domain[:, 1]) & all(pos .<= mesh.domain[:, 2])
#        # Get magnetic field at position by interpolation
#        fields, _ = grid(mesh, interpolator, pos)
#        bfield, _ = fields
#        # The movement should be dependent on the field direction, not the field
#        # strength 
#        B = norm(bfield) # field strength
#        b = bfield/B     # unit vector pointing i the magnetic field direction
#        @. pos = pos + b * stepsize
#        stepsforward += 1
#    end
#    # Backward
#    pos .= initpos
#    stepsbackward = -1 # We don't want to step if the the fist step is out of
#    # domain 
#    while all(pos .>= mesh.domain[:, 1]) & all(pos .<= mesh.domain[:, 2])
#        fields, _ = grid(mesh, interpolator, pos)
#        bfield, _ = fields
#        B = norm(bfield)
#        b = bfield/B
#        @. pos = pos - b * stepsize
#        stepsbackward += 1
#    end

    # Store line
    stepsforward = numsteps
    stepsbackward = numsteps
    lineforward  = zeros((mesh.numdims, stepsforward + 1))
    linebackward = zeros((mesh.numdims, stepsbackward + 1))
    lineforward[:, 1] .= initpos
    linebackward[:, 1] .= initpos

    # Follow field line forward
    for i = 1:stepsforward
        pos = lineforward[:, i]
        fields, _ = grid(mesh, interpolator, pos)
        bfield, _ = fields
        B = norm(bfield)
        b = bfield/B
        #@. lineforward[:, i + 1] = lineforward[:, i] + b * stepsize
        lineforward[:, i + 1] .= scheme(pos,
                                        stepsize,
                                        Solvers.fieldtracingforward,
                                        mesh.bField,
                                        interpolator,
                                        mesh.xCoords,
                                        mesh.yCoords,
                                        mesh.zCoords
                                        )
    end
    # Follow field line backward
    for i = 1:stepsbackward
        pos = linebackward[:, i]
        fields, _ = grid(mesh, interpolator, pos)
        bfield, _ = fields
        B = norm(bfield)
        b = bfield/B
        #@. linebackward[:, i + 1] = linebackward[:, i] - b * stepsize
        linebackward[:, i + 1] .= scheme(pos,
                                        stepsize,
                                        Solvers.fieldtracingbackward,
                                        mesh.bField,
                                        interpolator,
                                        mesh.xCoords,
                                        mesh.yCoords,
                                        mesh.zCoords
                                        )
    end

    # Fix points outside the domain
    for i = 1:stepsforward
        if any(lineforward[:,i] .< mesh.domain[:,1]) |
            any(lineforward[:,i] .> mesh.domain[:,2])
            lineforward[:, i:end] .= lineforward[:,i-1]
            break
        end
    end
    for i = 1:stepsbackward
        if any(linebackward[:,i] .< mesh.domain[:,1]) |
            any(linebackward[:,i] .> mesh.domain[:,2])
            linebackward[:, i:end] .= linebackward[:,i-1]
            break
        end
    end
    
    # Reverse linebackward and contatenate the lines
    fieldline = [linebackward[:, end:-1:1];; lineforward]
    return fieldline, lineforward, linebackward
end 


function plotfieldlines(mesh, numlines, stepsize, interpolator, scheme)
    fieldstrength = norm4(mesh.bField)
    maxfieldstrength = maximum(fieldstrength)
    # Create a function which returns the field strength at a given position.
    B(pos) = grid(fieldstrength,
                  interpolator,
                  pos,
                  mesh.xCoords,
                  mesh.yCoords,
                  mesh.zCoords)
    # Sample initial position for magnetic field lines based on the field
    # strength 
    positions, ratio = rejectionsampling(B,
                                         maxfieldstrength,
                                         numlines,
                                         mesh.domain)
    # Assumes the lines are no greater than twice the maximum domain extent.
    maxextent = maximum(mesh.domain[:, 2] .- mesh.domain[:, 1]) 
    numsteps = wpInt(2maxextent/stepsize)
    lines = zeros((mesh.numdims, numlines, 2*numsteps+2))
    for i = 1:numlines
        l, _, _ = generatefieldline(mesh,
                                    positions[:,i],
                                    stepsize,
                                    numsteps,
                                    interpolator,
                                    scheme
                                    )
        lines[:,i,:] .= l
    end
    PyPlot.figure()
    PyPlot.plot(lines[1,1,:], lines[2,1,:], color="black", linewidth=0.5)
    for i = 2:numlines
        PyPlot.plot(lines[1,i,:], lines[2,i,:], color="black", linewidth=0.5)
    end
    return lines, positions
end

end # module TPplots
