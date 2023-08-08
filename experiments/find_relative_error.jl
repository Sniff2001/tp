using TraceParticles
using Schemes
using Utilities
using PlasmoidParameters
using DipoleParameters
#using SpeiserParameters
using MirroringParameters
using Particles:                specieTable
using Meshes
using Interpolations_tp:        gridinterp, trilinear, trilinear_ip, quadratic_bspline, cubic_bspline
using LinearAlgebra

import Plots
using TPplots
import PyPlot

const plt = PyPlot
plt.pygui(true)


function runExperiment(
        experiment;                       # "Experiment"Parameters.jl function to run
        n_tuple_arr = [(1000,1000,2)],    # Grid size tuple (kept in array form for consistency)
        dt = [0.001],                     # Time step (same as above)
        run::Bool=false,                  # Whether or not to run the experiment or only make parameters
        save::Bool=true,                  # Whether or not to save the experiment (only works if run is also true)
        interpolator::String="trilinear", # Which interpolator to use
        laststep::Bool=true               # Last step or not
        )
        # Defining the arrays which contain the values of each experiment regarding position and velocity
        if run
                if laststep
                        pos_array = zeros(Float32, length(n_tuple_arr), length(dt), ndims, nparticles, 1)
                        vel_array = zeros(Float32, length(n_tuple_arr), length(dt), ndims, nparticles, 1)
                               
                else
                        pos_array = zeros(Float32, length(n_tuple_arr), length(dt), ndims, nparticles, nsteps+1)
                        vel_array = zeros(Float32, length(n_tuple_arr), length(dt), ndims, nparticles, nsteps+1)
                end
        end
        # Defining the array which contains the parameters
        parameter_array = Array{Parameters}(undef, length(n_tuple_arr), length(dt))
        for i in eachindex(n_tuple_arr)
                for j in eachindex(dt)
                        # Creating the parameters
                        params = experiment(n_tuple_arr[i], dt[j], interpolator=interpolator)
                        parameter_array[i,j] = params
                        if run
                                # Initialise experiment
                                exp = tp_init!(params)
                                # Run experiment
                                tp_run!(exp)
                                # Adds the current experiment to the arrays
                                if laststep
                                        pos_array[i,j,:,:] = exp.patch.tp.pos[:,:,end]
                                        vel_array[i,j,:,:] = exp.patch.tp.vel[:,:,end]
                                else
                                        pos_array[i,j,:,:,:] = exp.patch.tp.pos[:,:,:]
                                        vel_array[i,j,:,:,:] = exp.patch.tp.vel[:,:,:]
                                end
                                if save
                                        # Save simulation
                                        tp_save(exp)
                                end
                        end
                end # for loop
        end # for loop
        return pos_array, vel_array, parameter_array
end # function

function calculateError(
        standard,                               # Truth value to compare against
        nparticles::Integer;                    # Amount of particles in the experiment (note to self: can use parameters)
        pos_array=false,                        # Position array for all experiments
        vel_array=false,                        # Velocity array for all experiments
        parameter_array=false,                  # Parameters array for all experiments
        dt=false,                               # Time step array 
        error_type::String="relative error",    # Error type to use
        ndims::Integer=3,                       # Number of dimensions
        loadmode::Bool=false,                   # Whether or not to load the saved experiments
        plot_mode=plt.plot,                     # Plots the data for any valid PyPlot plotting
        show=true,
        laststep::Bool=true,                    # Last step or not
        label=""
        )
        error_dict = Dict(
                "absolute error" => (actual, measured) -> abs.(measured .- actual),
                "relative error" => (actual, measured) -> 100 .* abs.(measured ./ actual .- 1),
                # Below need a revamp of code before usable
                "chi squared" => (actual, measured) -> sum.((measured .- actual).^2 ./ actual),
                "root mean square" => (actual, measured) -> sqrt.(sum.((measured .- actual).^2) ./ nsteps)
        )
        if loadmode
                if laststep
                else
                        pos_array = zeros(Float32, length(n_tuple_arr), length(dt), ndims, nparticles, nsteps)
                        vel_array = zeros(Float32, length(n_tuple_arr), length(dt), ndims, nparticles, nsteps)
                
                end
                for i in axes(parameter_array, 1)
                        for j in axes(parameter_array, 2)
                                pos, vel = tp_loadtp(parameter_array[i,j], "$(parameter_array[i,j].tp_expdir)/$(parameter_array[i,j].tp_expname).tp")
                                pos_array[i,j,:,:] = pos[:,:,end]
                                vel_array[i,j,:,:] = vel[:,:,end]
                        end
                end
        end
        
        if length(n_tuple_arr) > 1
                for i in 1:nparticles
                        r_basis = mapslices(Utilities.norm, standard[:,1,:,i], dims=2)
                        r_norm = mapslices(Utilities.norm, pos_array[:,1,:,i], dims=2)
                        #display(standard)
                        #display(r_basis)
                        #println(size(r_basis), size(r_norm), size(standard), size(pos_array))
                        #error()
                        plt.subplot(2,2, i)
                        plot_mode(n, error_dict[error_type](r_basis, r_norm)[:], label="Particle $i ($label)")
                        if show
                                plt.xlabel("Grid size [N]")
                                plt.title("The variance in solution from grid size")
                                plt.ylabel(error_type)
                                plt.grid(ls=":")
                                plt.legend()
                                plt.show()
                        end
                end
        end
        if length(dt) > 1
                for i in 1:nparticles
                        r_basis = mapslices(Utilities.norm, standard[1,:,:,i], dims=2)
                        r_norm = mapslices(Utilities.norm, pos_array[1,:,:,i], dims=2)
                        #println(size(r_basis), size(r_norm), size(standard), size(pos_array))
                        plot_mode(dt, error_dict[error_type](r_basis, r_norm)[:], label="Particle $i ($label)")
                        plt.xlabel("dt [s]")
                        plt.title("The variance in solution from timestep")
                        plt.ylabel(error_type)
                        plt.grid(ls=":")
                        plt.legend()
                        plt.show()
                end
        end
end

function z(t, μ, mass, B0, L, A, ϕ)
        ω = √(2μ*B0/(mass*L^2))
        return @. A*sin(ω*t + ϕ)
end

function interpolationError(interpolator_results, labels; plot_mode=plt.plot)
        
        B0 = 30.0 # Magnetic field strenth parameter
        numparticles = 1
        species = 4*ones(Int32, numparticles) # Specifies the species of the particles
        L = 0.4
        a = 2.1L
        xi0 = (-a, -a, -a)
        # Upper bound of the three spatial axes
        xif = (a, a, a)
        xx, yy, zz, dx, dy, dz = Utilities.createaxes(xi0, xif, (100,100,100))
        numdims = 3
        Ex = 0.0
        Ey = 0.0
        Bfield = zeros(Float32, numdims, 100, 100, 100)
        Efield = zeros(size(Bfield))
        Efield[1,:,:,:] .= Ex
        Efield[2,:,:,:] .= Ey
        #Efield[1,:,1:30,:] .= 0.0
        Utilities.discretise!(Bfield, xx, yy, zz, Utilities.mirroringfield, B0, L)
        mass = specieTable[species[1], 1]
        charge = specieTable[species[1], 2]
        vel0 = [0.0, 0.1, 0.1]
        rL = mass*√(vel0[1]^2 + vel0[2]^2)/(charge*B0)
        pos0 = [-rL, 0.0, 0.0]
        mesh = Mesh(Bfield, Efield, xx, yy, zz)
        (B⃗, E⃗), _ = gridinterp(mesh,
                        trilinear,
                        pos0)
        B = norm(B⃗)
        b̂ = B⃗/B
        v = norm(vel0)
        vparal = [vel0 ⋅ b̂]
        vperp = √(v^2 - vparal[1]^2)
        μ = [mass*vperp^2/(2B)]
        Bmax = B0*norm(vel0)^2/(vel0[1]^2 + vel0[2]^2)
        zmax = L*√(Bmax/B0 - 1)
        ϕ = 0
        A = zmax
        time = 30
        times = collect(range(0.0, step=dt, length=round(Int, time/dt)+1))
        z_anal = z(times, μ[1], mass, B0, L, A, ϕ)

        i = 1
        for result in interpolator_results
                #println(size(n), size(result[:,1,3,:]), size(z_anal))
                #error()
                plot_mode(n, 100 .* abs.(result[:,1,3,:] ./ z_anal[end] .- 1), label=labels[i])
                i += 1
        end # for loop
        plt.xlabel("Grid size [N]")
        plt.title("The variance in solution from grid size")
        plt.ylabel("relative error")
        plt.grid(ls=":")
        plt.legend()
        plt.show()
end # interpolationError

function interpolatedFieldError(analytic_field, interpolated_field, particles, interpolation_label; plot_mode=plt.plot)
        interpolated_field_array = zeros(length(particles[1,1,1,1,:]))
        analytic_field_array = zeros(length(interpolated_field_array))
        for i in 1:length(particles[1,1,1,:,1])
                label = "Particle $i (" * interpolation_label * ")"
                plt.subplot(2,2, i)
                for t in 1:length(particles[1,1,1,1,:])
                        x = particles[1, 1, 1, i, t]
                        y = particles[1, 1, 2, i, t]
                        z = particles[1, 1, 3, i, t]
                        interpolated_field_array[t] = Utilities.norm(interpolated_field(x,y,z))
                        println("$label  $(round(100*t/length(particles[1,1,1,1,:]), 2))%")
                        analytic_field_array[t] = Utilities.norm(analytic_field(x,y,z))
                end
                plot_mode(collect(1:(length(particles[1,1,1,1,:]))), 100 .* abs.(interpolated_field_array ./ analytic_field_array .- 1), label=label)
                plt.xlabel("Step N")
                plt.title("The interpolation error compared to the analytical solution")
                plt.ylabel("relative error")
                plt.grid(ls=":")
                plt.legend()
        end
        #plt.show()
end

### If you want a constant time step or grid size
dt = 0.001                     # Plasmoid / Mirroring
#dt = 5.e-4                      # Dipole
#dt = 0.2e-8                     # Speiser
n_tuple_arr = [(100,100,2)]    # Plasmoid
#n_tuple_arr = [(256,256,256)]  # Dipole
#n_tuple_arr = [(100,100,100)]  # Speiser / Mirroring

### Some relevant variables
time = 2.9              # Simulation time
ndims=3                 # Dimensions of grid size
nparticles=4            # Number of particles

### Creates a logarithmic grid size array (assumes constant time step)
#n_full = unique(round.(Int, 10 .^ LinRange(0,3, 5)))
#n = n_full[n_full .> 1]                                 # Removes duplicate grid sizes
#n_tuple_arr = [(x,x,2) for x in n]                      # Turns it into the wanted format

### Creates a logarithmic time step array (assumes constant grid size)
#dt = unique(Float32, 10 .^ LinRange(-5,0, 60))


### Also used in the calculateError(), although I just defined it here for some reason
nsteps = round.(Int, time ./ dt)

### Creates an array of all the different positions and velocities from the changing variable
pos_array1, vel_array1, parameter_array1 = runExperiment(plasmoidParameters, n_tuple_arr=n_tuple_arr, dt=dt, run=true, save=false, interpolator="trilinear-ip", laststep=false)
pos_array2, vel_array2, parameter_array2 = runExperiment(plasmoidParameters, n_tuple_arr=n_tuple_arr, dt=dt, run=true, save=false, interpolator="quadratic-bspline-ip", laststep=false)
pos_array3, vel_array3, parameter_array3 = runExperiment(plasmoidParameters, n_tuple_arr=n_tuple_arr, dt=dt, run=true, save=false, interpolator="cubic-bspline-ip", laststep=false)
# last step = true for everything except interpolatedFieldError

#interpolationError([vel_array1, vel_array2, vel_array3], ["trilinear", "quadratic bspline", "cubic bspline"], plot_mode=plt.loglog)

### Creates a standard "truth" value to compare the error against
#standard, standard_vel, standard_parameter = runExperiment(plasmoidParameters, n_tuple_arr=[(1000,1000,2)], dt=[1e-3], run=true, save=false, interpolator="trilinear-ip")

### Some draft of getting standard out from the first call
#"""
#pos1, vel1 = tp_loadtp(standard_parameter[1,1], "$(standard_parameter[1,1].tp_expdir)/$(standard_parameter[1,1].tp_expname).tp")
#standard[:,:] = pos1[:,:,end]
#standard_vel[:,:] = vel1[:,:,end]
#"""

### Calculates and plots the desired error
#calculateError(standard, nparticles, pos_array=pos_array1, vel_array=vel_array1, parameter_array=parameter_array1, dt=dt, ndims=ndims, error_type="relative error", loadmode=false, plot_mode=plt.loglog, show=false, label="trilinear")
#calculateError(standard, nparticles, pos_array=pos_array2, vel_array=vel_array2, parameter_array=parameter_array2, dt=dt, ndims=ndims, error_type="relative error", loadmode=false, plot_mode=plt.loglog, show=false, label="quadratic bspline")
#calculateError(standard, nparticles, pos_array=pos_array3, vel_array=vel_array3, parameter_array=parameter_array3, dt=dt, ndims=ndims, error_type="relative error", loadmode=false, plot_mode=plt.loglog, show=true, label="cubic bspline")

### Calculates and plots error of field interpolation
xi0 = (0., 0., 0.)
xif = (1., 1., 1.)
bamp   = 10.0            # Amplification factor
bconst = [100., 0., 0]   # Constant term
bz     = 0.0             # z-component
μx = (xif[1] - xi0[1])/2 
μy = (xif[2] - xi0[2])/2 
σx = (xif[1] - xi0[1])/8
σy = (xif[2] - xi0[2])/8
μ  = (μx, μy)
σ  = (σx, σy)

xx, yy, zz, dx, dy, dz = Utilities.createaxes(xi0, xif, (100,100,2))
domainaxes, gridsizes, A = Utilities.normal3Donlyz(xi0, xif, (100,100,2), μ, σ)
Bfield_interp = Schemes.curl(A, gridsizes, Schemes.derivateCentral)
@. Bfield_interp = bamp*Bfield_interp + bconst
@. Bfield_interp[3,:,:,:] = bz
Ex = 0.0
Ey = 0.0
Ez = 50.0
Efield_interp = zeros(Float64, size(Bfield_interp))
Efield_interp[1, :, :, :] .= Ex
Efield_interp[2, :, :, :] .= Ey
Efield_interp[3, :, :, :] .= Ez

trilinear_interp_fields = trilinear_ip(Bfield_interp, Efield_interp, xx, yy, zz, (false, false, true))
qbspline_interp_fields = quadratic_bspline(Bfield_interp, Efield_interp, xx, yy, zz, (false, false, true))
cbspline_interp_fields = cubic_bspline(Bfield_interp, Efield_interp, xx, yy, zz, (false, false, true))

analytic_Bfield = (x,y,z) -> bamp*Utilities.normal3Donlyz(μ, σ)[1](x,y) + bconst

interpolatedFieldError(analytic_Bfield, trilinear_interp_fields[1], pos_array1, "trilinear")
interpolatedFieldError(analytic_Bfield, qbspline_interp_fields[1], pos_array2, "quadratic bspline")
interpolatedFieldError(analytic_Bfield, cbspline_interp_fields[1], pos_array3, "cubic bspline")
plt.show()