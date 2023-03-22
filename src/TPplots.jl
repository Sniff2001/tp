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
const plt = PyPlot
using LinearAlgebra:    norm
using LaTeXStrings
# Internal libraries
using WorkingPrecision: wpInt, wpFloat
using Meshes
using Patches:          Patch
using Particles
using Utilities:        rejectionsampling, norm4, norm2
using Interpolations
using Schemes
using Solvers

export plt
export plot

#------#
# Mesh #
#-------------------------------------------------------------------------------
function quiverslice!(
    ax  ::plt.PyCall.PyObject,
    mesh::Mesh, 
    normal::String, 
    point::wpInt,
    field::String="B"
    )
    if field == "B"
        f = mesh.bField
    elseif field == "E"
        f = mesh.eField
    else
        println("Error: Invalid field.")
    end
    if normal == "x"
        ax.quiver(mesh.yCoords,
               mesh.zCoords,
               transpose(f[2, point, :, :]),
               transpose(f[3, point, :, :])
               )
    elseif normal == "y"
        ax.quiver(mesh.xCoords,
               mesh.zCoords,
               transpose(f[1, :, point, :]),
               transpose(f[3, :, point, :])  
               )
    elseif normal == "z"
        ax.quiver(mesh.xCoords,
               mesh.yCoords,
               transpose(f[1, :, :, point]),
               transpose(f[2, :, :, point])  
               )
    else
        println("Error: Plane not valid.")
    end
end # function quiverslice


function streamplotslice!(
    ax    ::plt.PyCall.PyObject,
    mesh  ::Mesh, 
    normal::String, 
    point ::wpInt,
    field ::String="B")
    if field == "B"
        f = mesh.bField
    elseif field == "E"
        f = mesh.eField
    else
        println("Error: Invalid field.")
    end
    if normal == "x"
        xx = mesh.yCoords
        yy = mesh.zCoords
        uu = transpose(f[2, point, :, :])
        vv = transpose(f[3, point, :, :])
    elseif normal == "y"
        xx = mesh.xCoords
        yy = mesh.zCoords
        uu = transpose(f[1, :, point, :])
        vv = transpose(f[3, :, point, :])
    elseif normal == "z"
        xx = mesh.xCoords
        yy = mesh.yCoords
        uu = transpose(f[1, :, :, point])
        vv = transpose(f[2, :, :, point])
    else
        println("Error: Plane not valid.")
    end
    ax.streamplot(
        xx,
        yy,
        uu,
        vv,
        linewidth=0.3,
        arrowsize=0.6,
        color="black"
    )
end # function streamplotslice


function pcolormeshslice!(
    ax    ::plt.PyCall.PyObject,
    mesh  ::Mesh, 
    normal::String, 
    point ::wpInt,
    field ::String="B")
    if field == "B"
        f = mesh.bField
        label = latexstring("\$|\\mathbf{B}|\$")
    elseif field == "E"
        f = mesh.eField
        label = latexstring("\$|\\mathbf{E}|\$")
    else
        println("Error: Invalid field.")
    end
    fabs = norm4(f)
    if normal == "x"
        pfabs = copy(transpose(fabs[point,:,:]))
        uuu = mesh.yCoords' .* ones(length(mesh.zCoords))
        vvv = ones(length(mesh.yCoords))' .* mesh.zCoords
    elseif normal == "y"
        pfabs = copy(transpose(fabs[:,point,:]))
        uuu = mesh.xCoords' .* ones(length(mesh.zCoords))
        vvv = ones(length(mesh.xCoords))' .* mesh.zCoords
    elseif normal == "z"
        pfabs = copy(transpose(fabs[:,:,point]))
        uuu = mesh.xCoords' .* ones(length(mesh.yCoords))
        vvv = ones(length(mesh.xCoords))' .* mesh.yCoords
    else
        println("Error: Plane not valid.")
    end
    pcm = ax.pcolormesh(uuu, vvv, pfabs,
                        alpha=1.0,
                        cmap=plt.get_cmap("Greys"))
    cb = plt.colorbar(pcm,
                      label=label,
                      ax=ax)
end # function colormeshslice


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



function plotperiodictrajectory(
    ax::plt.PyCall.PyObject,
    pos::Matrix{wpFloat},
    partidx::wpInt,
    domain ::Matrix{wpFloat},
    cm,
    colorrange,
    )
    extent = domain[:,2] .- domain[:,1]
    println(extent)
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
            ax.plot(pos[1,j:i], pos[2,j:i], color=cm(colorrange[partidx]),
                        label="Particle $partidx")                       
            j = i+1
            break
        end
        ax.plot(pos[1,j:end], pos[2,j:end], color=cm(colorrange[partidx]))
    else
        ax.plot(pos[1,:], pos[2,:], color=cm(colorrange[partidx]),
                    label="Particle $partidx")
    end
end 


#---------#
# Patches #
#-------------------------------------------------------------------------------
function plot(
    patch    ::Patch,
    labelling::Bool=false
    )

    pos = getpos(patch.tp)
    times = collect(range(0.0, step=patch.dt, length=patch.numSteps+1))
    Ek = kineticenergy(patch.tp)

    fig, axes = plt.subplots(1,1)

    # Make streamplot of magnetic field
    streamplotslice!(axes, patch.mesh, "z", 1)
    # Make pcolormesh of magnetic field strength
    pcolormeshslice!(axes, patch.mesh, "z", 1)
    # Plot particle trajectories
    cm = plt.get_cmap(:tab20)
    colorrange = (0:patch.numParticles) ./ patch.numParticles
    for i = 1:patch.numParticles
        plotperiodictrajectory(axes, pos[:,i,:], i, patch.mesh.domain,
                               cm, colorrange)
    end
    # Mark initial positions of particles
    for i = 1:patch.numParticles
        if i == 1
            axes.plot(pos[1,1,1], pos[2,1,1], marker=".", color="Black",
                        label=latexstring("\$t_0\$"), linestyle="None")
        else
            axes.plot(pos[1,i,1], pos[2,i,1], marker=".", color="Black")
        end
    end

    # Set title and labelling
    if typeof(patch.tp) == ParticleSoA
        axes.set_title("Full-orbit")
    else
        axes.set_title("GCA")
    end
    if labelling
        axes.legend()
    end
    axes.set_xlabel(latexstring("\$x\$"))
    axes.set_ylabel(latexstring("\$y\$"))
   
    # Make energy distribution-plots
    numbins = 20
    p1 = plotenergydistr(patch.tp, 1, numbins, "Initial energy")
    p2 = plotenergydistr(patch.tp, patch.numSteps+1, numbins, "final energy")
    Plots.plot(p1,p2)#p3,p4)

end # function plot



#---------------------#
# Field line tracing  #
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


function plotfieldlines(mesh, numlines::wpInt, stepsize, interpolator, scheme)
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
    plt.figure()
    plt.plot(lines[1,1,:], lines[2,1,:], color="black", linewidth=0.5)
    for i = 2:numlines
        plt.plot(lines[1,i,:], lines[2,i,:], color="black", linewidth=0.5)
    end
    return lines, positions
end
#|
function plotfieldlines(
    mesh,
    initpos::Matrix{wpFloat},
    stepsize,
    interpolator,
    scheme
    )
    numlines = size(initpos)[2]
    positions = initpos
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
    plt.figure()
    plt.plot(lines[1,1,:], lines[2,1,:], color="black", linewidth=0.5)
    for i = 2:numlines
        plt.plot(lines[1,i,:], lines[2,i,:], color="black", linewidth=0.5)
    end
    return lines, positions
end


#---------------------------------------------------#
# Messy plotting function hardcoded for RoCMI 2023  #
#-------------------------------------------------------------------------------
"""
    plotRoCMI(patch, labelling)
Plotting function used to make RoCMI-poster figures.
"""
function plotRoCMI(
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


"""
    plotRoCMI(patch, labelling)
Function for attempt to make nice subplots for the RoCMI-poster (Did not finish).
"""
function plotRoCMIsubplot(
    fbpatch  ::Patch,
    GCApatch ::Patch,
    labelling::Bool=false
    )
    times = collect(range(0.0, step=fbpatch.dt, length=fbpatch.numSteps+1))
    Ek = kineticenergy(fbpatch.tp)
    plotplasmoid(
        fbpatch.tp,
        GCApatch.tp,
        times,
        Ek,
        fbpatch.mesh.xCoords, fbpatch.mesh.yCoords, fbpatch.mesh.zCoords,
        fbpatch.mesh.bField,
        fbpatch.numParticles,
        labelling
    )
end


"""
    plotplasmoid(tp::TraceParticle, ...)
Function for plotting GCA or full orbit trajectories of the plasmoid-experiment
used to make figures for the RoCMI 2023 poster.
"""
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
    plt.figure()
    for i = 1:numparticles
        if i == 2
            plt.plot(times[1:1904], Ek[i, 1:1904], label="Particle $i")
        elseif i == 3
            plt.plot(times[1:1744], Ek[i, 1:1744], label="Particle $i")
        elseif i == 4
            plt.plot(times[1:793], Ek[i, 1:793], label="Particle $i")
        else
            plt.plot(times, Ek[i, :], label="Particle $i")
        end
    end
    if typeof(tp) == ParticleSoA
        plt.title("Full-orbit: Kinetic energy")
    else
        plt.title("GCA: kinetic energy")
    end
    plt.xlabel("Time, s")
    plt.ylabel("Energy, J")
    if labelling
        plt.legend()
    end
    
    # Streamplot of magnetic field with initial positions and velocity
    #plt.figure()
    #plt.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]))
    #for i = 1:numparticles
    #    plt.plot(pos[1,i,1], pos[2,i,1], marker="o")
    #end
    #if typeof(tp) == ParticleSoA
    #    # Plot starting velocity as an arrow
    #    vel = getvel(tp)
    #    plt.quiver(pos[1, :, 1], pos[2, :, 1],
    #                  vel[1, :, 1], vel[2, :, 1],
    #                  width=0.003)
    #    plt.title("Full-orbit trajectories")
    #else
    #    plt.title("GCA trajectories")
    #end
    
    # Streamplot of magnetic field with particle trajectories
    babs = norm4(B) # Magnetic field strength
    plt.figure()
    #plt.axis("equal")
    #plt.contourf(xx, yy, transpose(babs[:,:,1]), levels=70,
    #                alpha=1.0,
    #                cmap=plt.get_cmap("Greys"))
    pbabs = copy(transpose(babs[:,:,1]))
    uuu = xx' .* ones(length(yy))
    vvv = ones(length(xx))' .* yy
    plt.pcolormesh(uuu, vvv, pbabs,
                    alpha=1.0,
                    cmap=plt.get_cmap("Greys"))
    cb = plt.colorbar(label=latexstring("\$|\\mathbf{B}|\$"),
                         ticks=[minimum(pbabs), maximum(pbabs)]
                         )
    println(minimum(pbabs))
    cb.ax.set_yticklabels(
        [latexstring("\$B_{min}\$"),latexstring("\$B_{max}\$")])
    plt.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]),
                      linewidth=0.3,
                      arrowsize=0.6,
                      color="black")
    plt.xlabel(latexstring("\$x\$"))
    plt.ylabel(latexstring("\$y\$"))
    cm = plt.get_cmap(:tab20)
    colorrange = (0:numparticles) ./ numparticles

    for i = 1:numparticles
        plotperiodictrajectory(pos[:,i,:], i, cm, colorrange)
    end
    for i = 1:numparticles
        # Mark position after t=1.5
        if i == 1
            #plt.plot(pos[1,1,1501], pos[2,1,1501], marker="^",
            #            color="black",label="t = 1.5s", linestyle="None") 
            plt.plot(pos[1,1,1], pos[2,1,1], marker=".", color="Black",
                        label=latexstring("\$t_0\$"), linestyle="None")
        else
            #plt.plot(pos[1,i,1501], pos[2,i,1501], marker="^",
            #            color="black") 
            plt.plot(pos[1,i,1], pos[2,i,1], marker=".", color="Black")
        end
    end

    if typeof(tp) == ParticleSoA
        # Plot starting velocity as an arrow
        #vel = getvel(tp)
        #plt.quiver(pos[1, :, 1], pos[2, :, 1],
        #              vel[1, :, 1], vel[2, :, 1],
        #              width=0.003)
        plt.title("Full-orbit")
    else
        plt.title("GCA")
    end

    if labelling
        plt.legend()
    end

end # function plotplasmoid


"""
    plotplasmoid(tp::ParticleSoA, tp::GCAParticleSoA, ...)
Subplot both full orbit and GCA for comparison in RoCMI-poster.
"""
function plotplasmoid(
    fbtp::ParticleSoA,
    gcatp::GCAParticleSoA,
    times,
    Ek,
    xx,
    yy, 
    zz,
    B,
    numparticles,
    labelling
    )
    #fig = plt.figure(figsize=(10,6))
    fig = plt.figure()
    energyplot = fig.add_subplot(3,4,(9,12))
    # Energy plot
    for i = 1:numparticles
        if i == 2
            energyplot.plot(times[1:1904], Ek[i, 1:1904], label="Particle $i")
        elseif i == 3
            energyplot.plot(times[1:1744], Ek[i, 1:1744], label="Particle $i")
        elseif i == 4
            energyplot.plot(times[1:793], Ek[i, 1:793], label="Particle $i")
        else
            energyplot.plot(times, Ek[i, :], label="Particle $i")
        end
    end
    energyplot.set_title("Full-orbit: Kinetic energy")
    energyplot.set_xlabel("Time, s")
    energyplot.set_ylabel("Energy, J")
    if labelling
        energyplot.legend()
    end
    
    # Streamplot of magnetic field with particle trajectories
    fbplot = fig.add_subplot(3,4,(3))
    pos = fbplot.get_position()
    fbplot.set_position([pos.x0*2, pos.y0*2, pos.width, pos.height/2])
    gcaplot = fig.add_subplot(3,4,(1))
    pos = gcaplot.get_position()
    fbplot.set_title("Full-orbit")
    gcaplot.set_title("GCA")
    fbplot.set_xlabel(latexstring("\$x\$"))
    fbplot.set_ylabel(latexstring("\$y\$"))
    gcaplot.set_xlabel(latexstring("\$x\$"))
    # Full orbit
    babs = norm4(B) # Magnetic field strength
    pcm1 = fbplot.contourf(xx, yy, transpose(babs[:,:,1]),
                    alpha=1.0,
                    cmap=plt.get_cmap("Greys"))
    fbcb = plt.colorbar(pcm1,
#                         label=latexstring("\$|\\mathbf{B}|\$"),
                         ticks=[0, maximum(babs)],
                         ax=fbplot
                         )
    fbcb.ax.set_yticklabels(["0",latexstring("\$B_{max}\$")])
    fbplot.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]),
                      linewidth=0.3,
                      arrowsize=0.6,
                      color="black")
    # GCA
    pcm2 = gcaplot.contourf(xx, yy, transpose(babs[:,:,1]),
                           alpha=1.0,
                           cmap=plt.get_cmap("Greys"))
    
    gcacb = plt.colorbar(pcm2,
#                         label=latexstring("\$|\\mathbf{B}|\$"),
                         ticks=[0, maximum(babs)],
                         ax=gcaplot
                         )
    gcacb.ax.set_yticklabels(["0",latexstring("\$B_{max}\$")])
    gcaplot.streamplot(xx, yy, transpose(B[1,:,:,1]), transpose(B[2,:,:,1]),
                      linewidth=0.3,
                      arrowsize=0.6,
                      color="black")

    # Plot particle trajectories
    cm = plt.get_cmap(:tab20)
    colorrange = (0:numparticles) ./ numparticles
    for i = 1:numparticles
        plotperiodictrajectory(fbplot, fbtp.pos[:,i,:], i, cm, colorrange)
        plotperiodictrajectory(gcaplot, gcatp.R[:,i,:], i, cm, colorrange)
    end

    for i = 1:numparticles
        # Mark position after t=1.5
        if i == 1
            fbplot.plot(fbtp.pos[1,1,1], fbtp.pos[2,1,1], marker=".", color="Black",
                        label=latexstring("\$t_0\$"), linestyle="None")
            gcaplot.plot(gcatp.R[1,1,1], gcatp.R[2,1,1], marker=".", color="Black",
                        label=latexstring("\$t_0\$"), linestyle="None")
        else
            fbplot.plot(fbtp.pos[1,i,1], fbtp.pos[2,i,1], marker=".", color="Black")
            gcaplot.plot(gcatp.R[1,i,1], gcatp.R[2,i,1], marker=".", color="Black")
        end
    end

    if labelling
        fbplot.legend()
        gcaplot.legend()
    end

end # function plotplasmoid


function plotperiodictrajectoryRoCMI(pos, partidx, cm, colorrange)
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
            plt.plot(pos[1,j:i], pos[2,j:i], #color=cm(colorrange[partidx]),
                        label="Particle $partidx")                       
            j = i+1
            break
        end
    else
        plt.plot(pos[1,:], pos[2,:], #color=cm(colorrange[partidx]),
                    label="Particle $partidx")
    end
end 
#|
function plotperiodictrajectoryRoCMI(pyaxes, pos, partidx, cm, colorrange)
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
            pyaxes.plot(pos[1,j:i], pos[2,j:i], #color=cm(colorrange[partidx]),
                        label="Particle $partidx")                       
            j = i+1
            break
        end
    else
        pyaxes.plot(pos[1,:], pos[2,:], #color=cm(colorrange[partidx]),
                    label="Particle $partidx")
    end
end 



end # module TPplots
