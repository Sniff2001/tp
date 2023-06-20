#-------------------------------------------------------------------------------
# Created 05.06.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                TraceParticles.jl
#
#-------------------------------------------------------------------------------
# The tp (trace particles) package.
#-------------------------------------------------------------------------------

module TraceParticles

using Mmap

using Bifrost

using Patches
using Meshes
using Particles
using Solvers
using Schemes
using Interpolations_tp
using Utilities
using Constants

export Parameters
export Experiment
export tp_init!
export tp_run!
export tp_save
export tp_load
export tp_loadtp
export tp_loadmesh
export tp_loadbg


methodmap = Dict(
    "full-orbit"     => Solvers.fullOrbit,
    "full-orbit-isf" => Solvers.fullOrbit_interstaticfield,
    "full-orbit-rel" => Solvers.relFullOrbitExplLeapFrog,
    "GCA"            => Solvers.GCA,
    "RK4"            => Schemes.rk4,
    "trilinear"      => Interpolations_tp.trilinear,
    "bilinear_xz"    => Interpolations_tp.bilinear_xz,
)


mutable struct Parameters
    #
    npart ::Integer
    nsteps::Integer
    dt    ::Real
    # Working precision on snap (mesh and fields) and part (trace particles)
    wp_snap::DataType
    wp_part::DataType
    # Solver, numerical scheme and interpolation scheme
    solver::String
    scheme::String
    interp::String
    # Location for saving simulation-data
    tp_expname::String
    tp_expdir ::String
    # Bifrost-snapshot
    br_expname::String
    br_expdir ::String
    br_isnap  ::Integer
    # Manual mesh and EM-fields
    x ::Vector{T} where {T<:Real}
    y ::Vector{T} where {T<:Real}
    z ::Vector{T} where {T<:Real}
    bx::Array{T, 3} where {T<:Real}
    by::Array{T, 3} where {T<:Real}
    bz::Array{T, 3} where {T<:Real}
    ex::Array{T, 3} where {T<:Real}
    ey::Array{T, 3} where {T<:Real}
    ez::Array{T, 3} where {T<:Real}
    # Number of grid points in the x, y, and z-direction
    nx::Integer
    ny::Integer
    nz::Integer
    # Particle-initialisation
    specie   ::Vector{T} where {T<:Integer}
    pos_distr::String # Particle position distribution "uniform": uniform
    vel_distr::String # Particle velocity distribution    ↑ or "mb": MB
    posxbounds::Vector{T} where {T<:Real} 
    posybounds::Vector{T} where {T<:Real} 
    poszbounds::Vector{T} where {T<:Real} 
    velxbounds::Vector{T} where {T<:Real} 
    velybounds::Vector{T} where {T<:Real} 
    velzbounds::Vector{T} where {T<:Real} 
    T         ::Real
    seed      ::Integer
    # Boundary conditions
    pbc::Tuple{Bool, Bool, Bool} # (x,y,z) If not periodic boundary conditions,
                                 # particles are
                                 # killed when crossing all boundaries.
    # "Private" parameters
    bg_input::String # "br": Bifrost snapshot
                     # "manual": Mesh and fields are defined manually in Julia.

    # Constructors
    #--------------------------------------------------------------------------
    # Bifrost-snap
    function Parameters(
        ;
        # Required parameters
        #-----------------------------------------|
        npart ::Integer,
        nsteps::Integer,
        dt    ::Real,
        # Working precision
        wp_snap::DataType,
        wp_part::DataType,
        #
        # Optionals
        #-----------------------------------------|
        solver::String="full-orbit",
        scheme::String="RK4",
        interp::String="trilinear",
        tp_expname=nothing,
        tp_expdir=nothing,
        # Either Bifrost-snapshot
        br_expname=nothing,
        br_expdir =nothing,
        br_isnap  =nothing,
        # or manual mesh and fields
        x =nothing,
        y =nothing,
        z =nothing,
        bx=nothing,
        by=nothing,
        bz=nothing,
        ex=nothing,
        ey=nothing,
        ez=nothing,
        # Either way, gridsize will be useful
        nx=nothing,
        ny=nothing,
        nz=nothing,
        # Particle-initialisation
        specie=ones(wp_part, npart),
        pos_distr::String="uniform",
        vel_distr::String="point",
        # Default value of bounds depend on the mesh which may be given by
        # a Bifrost snapshot.
        posxbounds=nothing,
        posybounds=nothing,
        poszbounds=nothing,
        velxbounds=zeros(wp_part, 2),
        velybounds=zeros(wp_part, 2),
        velzbounds=zeros(wp_part, 2),
        T         =nothing,
        seed      =wp_part(0),
        # Periodic boundary conditions
        pbc::Tuple{Bool, Bool, Bool}=(false, false, false),
        # "Private" parameters
        #-----------------------------------------|
        bg_input::String="br"
        ) 
        #
        println("tp.jl: Constructing Parameters...")
        #
        #  Create parameter instance with all parameters undefined
        params = new()

        #-----------------------------------------------------------------------
        # Check if mesh and EM-field paramters are present.
        # Bifrost snapshot is preferred by default.
        if all(i -> i != nothing, (br_expname, br_expdir, br_isnap))
            params.bg_input = "br"
            params.br_expname = br_expname
            params.br_expdir = br_expdir
            params.br_isnap = br_isnap
        elseif all(i -> i != nothing, (x, y, z, bx, by, bz, ex, ey, ez))
            params.bg_input = "manual"
            params.x = x
            params.y = y
            params.z = z
            params.nx = Int32(length(x))
            params.ny = Int32(length(y))
            params.nz = Int32(length(z))
            params.bx = bx
            params.by = by
            params.bz = bz
            params.ey = ex
            params.ex = ey
            params.ez = ez
        else
            error("Missing parameter(s): mesh or/and EM-field parameters")
        end 

        #-----------------------------------------------------------------------
        # Check if optional paramters without default are available
        if tp_expname != nothing
            params.tp_expname = tp_expname
        end
        if tp_expdir != nothing
            params.tp_expdir = tp_expdir
        end
        if posxbounds != nothing
            params.posxbounds = posxbounds
        end
        if posybounds != nothing
            params.posybounds = posybounds
        end
        if poszbounds != nothing
            params.poszbounds = poszbounds
        end
        if nx != nothing
            params.nx = nx
        end
        if ny != nothing
            params.ny = ny
        end
        if nz != nothing
            params.nz = nz
        end
        if T != nothing
            params.T = T
        end

        #-----------------------------------------------------------------------
        # Set required parameters
        params.npart = npart
        params.nsteps = nsteps
        params.dt = dt
        params.wp_snap = wp_snap
        params.wp_part = wp_part
        # Set optional parameters
        params.solver = solver
        params.scheme = scheme
        params.interp = interp
        params.specie = specie
        params.pos_distr = pos_distr
        params.vel_distr = vel_distr
        params.seed = seed
        params.pbc = pbc
        params.velxbounds = velxbounds
        params.velybounds = velybounds
        params.velzbounds = velzbounds
        # Set default value to optional parameters that are not passed and not
        # already defined above.

        # Construct methodmap

        # Return Parameters instance
        return params
    end # Constructor Parameters 
end # struct Parameters


struct Experiment
    params::Parameters
    patch ::Patch
end


function tp_checkrequirements(
    params::Parameters
    )
    if params.npart == nothing
        error("Missing parameter: npart")
    elseif params.nsteps == nothing
        error("Missing parameter: nsteps")
    elseif params.dt == nothing
        error("Missing parameter: dt")
    elseif params.solver == nothing
        error("Missing parameter: solver")
    elseif params.scheme == nothing
        error("Missing parameter: scheme")
    elseif params.interp == nothing
        error("Missing parameter: interp")
    end
    if all(i -> isdefined(params, i), (:br_expname,
                                       :br_expdir,
                                       :br_isnap
                                       )
           )
    elseif all(i -> isdefined(params, i), (:x, :y, :z,
                                           :bx, :by, :bz,
                                           :ex, :ey, :ez
                                           )
               )
    else
        error("Missing parameter(s): mesh or/and EM-field")
    end
    if params.bg_input == nothing
        error("Missing paramter: bg_input")
    end
    if params.npart > length(params.specie)
        error("Number of particles outnumbers the length of specie-parameter.")
    end

end # function checkRequirements
    

function tp_init!(
    params::Parameters
    )
    # Check whether the required parameters are present. Since the struct is
    # mutable, it may have changed since constrcution.
    tp_checkrequirements(params)

    # Define variables to avoid magic numbers
    numdims = 3

    println("tp.jl: Initialising experiment...")

function tp_initfromfile(
    params::Parameters,
    mesh_filename::String,
    bg_filename::String,
    tp_filename::String,
    )
    #
    println("tp.jl: Initialising experiment from file...")
    #
    # Open tp-file
    #
    pos, vel = tp_loadtp(params, tp_filename)
    #
    # Open mesh-file
    #
    x, y, z = tp_loadmesh(params, mesh_filename)
    #
    # Open bg-file
    #
    bField, eField, ∇B, ∇b̂, ∇ExB = tp_loadbg(params, bg_filename) 
    #
    # Create Mesh, ParticleSoA, Patch and Experiment
    #
    mesh = Mesh(bField, eField, ∇B, ∇b̂, ∇ExB, x, y, z)
    particles = ParticleSoA(pos, vel, params.specie)
    patch = Patch(mesh,
                  particles,
                  methodmap[params.solver],
                  methodmap[params.scheme],
                  methodmap[params.interp],
                  params.dt,
                  params.nsteps,
                  params.npart,
                  params.pbc
                  )
    exp = Experiment(params, patch);
    #
    return exp, params
end # function tp_init
#|


function tp_loadtp(
    params::Parameters,
    filename::String
    )
    numdims = 3
    f = open(filename)
    pos = zeros(params.wp_part, numdims, params.npart, params.nsteps)
    vel = zeros(params.wp_part, numdims, params.npart, params.nsteps)
    pos[1,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    pos[2,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    pos[3,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    vel[1,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    vel[2,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    vel[3,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    close(f)
    return pos, vel
end


function tp_loadmesh(
    params::Parameters,
    filename::String
    )
    f = open(filename)
    x = mmap(f, Vector{params.wp_snap}, params.nx)
    y = mmap(f, Vector{params.wp_snap}, params.ny)
    z = mmap(f, Vector{params.wp_snap}, params.nz)
    close(f)
    return x, y, z
end
    

function tp_loadbg(
    params::Parameters,
    filename::String
    )
    meshsize = (params.nx, params.ny, params.nz)
    bField = zeros(params.wp_snap, 3, meshsize...)
    eField = zeros(params.wp_snap, 3, meshsize...)
    ∇B = zeros(params.wp_snap, 3, meshsize...)
    ∇b̂ = zeros(params.wp_snap, 3, 3, meshsize...)
    ∇ExB = zeros(params.wp_snap, 3, 3, meshsize...)
    f = open(filename)
    for i = 1:3
        bField[i,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
    end
    for i = 1:3
        eField[i,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
    end
    for i = 1:3
        ∇B[i,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
    end
    for i = 1:3
        for j = 1:3
            ∇b̂[i,j,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
        end
    end
    for i = 1:3
        for j = 1:3
            ∇ExB[i,j,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
        end
    end
    close(f)
    return bField, eField, ∇B, ∇b̂, ∇ExB
end
    #---------------------------------------------------------------------------
    # Construct mesh
    #---------------------------------------------------------------------------
    if params.bg_input == "br"
        mesh = Mesh(params.br_expname, params.br_expdir, params.br_isnap)
        params.nx = length(mesh.xCoords)
        params.ny = length(mesh.yCoords)
        params.nz = length(mesh.zCoords)
    else params.bg.input == "manual"
        meshsize = (params.nx, params.ny, params.nz)
        bField = zeros(wp, 3, meshsize...)
        eField = zeros(wp, 3, meshsize...)
        bField[1,:,:,:] = bx
        bField[2,:,:,:] = by
        bField[3,:,:,:] = bz
        eField[1,:,:,:] = ex
        eField[2,:,:,:] = ey
        eField[3,:,:,:] = ez
        mesh = Mesh(bField, eField, params.x, params.y, params.z)
    end 
    
    #---------------------------------------------------------------------------
    # Construct particles
    #---------------------------------------------------------------------------
    #
    # First, create positions and velocities according to given parameters
    #
    # Go thorugh the various intial distributions
    # Positions
    if params.pos_distr == "uniform"
        # Set default bounds if not defined
        # Default bounds are the mesh domain boundaries.
        if !isdefined(params, :posxbounds)
            params.posxbounds = [mesh.xCoords[1], mesh.xCoords[end]]
        end
        if !isdefined(params, :posybounds)
            params.posybounds = [mesh.yCoords[1], mesh.yCoords[end]]
        end
        if !isdefined(params, :poszbounds)
            params.poszbounds = [mesh.zCoords[1], mesh.zCoords[end]]
        end
        pos = Utilities.inituniform(params.npart,
                                    params.posxbounds,
                                    params.posybounds,
                                    params.poszbounds,
                                    params.wp_part,
                                    params.seed
                                    )
    elseif params.pos_distr == "point"
        if !isdefined(params, :posxbounds)
            params.posxbounds = [(mesh.xCoords[1] - mesh.xCoords[end])/2.0]
        end
        if !isdefined(params, :posybounds)
            params.posybounds = [(mesh.yCoords[1] - mesh.yCoords[end])/2.0]
        end
        if !isdefined(params, :poszbounds)
            params.poszbounds = [(mesh.zCoords[1] - mesh.zCoords[end])/2.0]
        end
        pos = np.ones(params.wp_part, numdims, params.npart)
        # If particles are given by coordinates
        if (length(params.posxbounds) ==
            length(params.posybounds) ==
            length(params.poszbounds) == params.npart)
            pos[1,:] = params.posxbounds
            pos[2,:] = params.posybounds
            pos[3,:] = params.poszbounds
        else
            pos[1,:] *= params.posxbounds[1]
            pos[2,:] *= params.posybounds[1]
            pos[3,:] *= params.poszbounds[1]
        end
    else
        error("Invalid parameter value: pos_distr")
    end
    # IMPLEMENT e.g
    #
    # position distributed according to density
    #

    # Velocities
    if params.vel_distr == "mb"
        # Get temperature 
        if params.bg_input == "br"
            vel = ones(numdims, params.npart)
            br_temp = dropdims(br_load_auxvariable(params.br_expname,
                                          [params.br_isnap],
                                          params.br_expdir,
                                          "tg",
                                          params.wp_snap
                                                   ),
                               dims=numdims + 1
                               )
            for i = 1:params.npart
                # Interpolate the temperature at the position of the particle
                x⃗ = pos[:,i]
                t, _ = gridinterp(br_temp,
                               methodmap[params.interp],
                               x⃗,
                               mesh.xCoords,
                               mesh.yCoords,
                               mesh.zCoords,
                               )

                # Standard deviation of velocity
                σ = √(Constants.k_B*t/specieTable[params.specie[i], 1])
                # Expectance value of particle velocity components is zero
                μ = 0.0
                # Draw velocity components from normal-distribution
                vel[:,i] = randn(μ, σ, (numdims))
            end
        else
            error(string("\"mb\" velocity distribution only available with ",
                         "Bifrost background. Use \"mb-onetemp\" to use one ",
                         "temperature for all particles"
                         )
                  )
        end
    elseif params.vel_distr == "mb-onetemp"
        if !isdefined(params, :T)
            error("Parameter not defined: T")
        end
        vel = ones(numdims, params.npart)
        for i = 1:params.npart
            # Standard deviation of velocity
            σ = √(Constants.k_B*params.T/specieTable[params.specie[i], 1])
            # Expectance value of particle velocity components is zero
            μ = 0.0
            # Draw velocity components from normal-distribution
            vel[:,i] = randn(μ, σ, (numdims))
        end
    elseif params.vel_distr == "uniform"
        vel = Utilities.inituniform(params.npart,
                                    params.velxbounds,
                                    params.velybounds,
                                    params.velzbounds,
                                    params.wp_part,
                                    params.seed
                                    )
    elseif params.vel_distr == "point"
        vel = ones(params.wp_part, numdims, params.npart)
        if (length(params.velxbounds) ==
            length(params.velybounds) ==
            length(params.velzbounds) == params.npart)
            vel[1,:] = params.velxbounds
            vel[2,:] = params.velybounds
            vel[3,:] = params.velzbounds
        else
            vel = ones(numdims, params.npart)
            vel[2,:] *= params.velybounds[1]
            vel[3,:] *= params.velzbounds[1]
        end
    else
        error("Invalid parameter value: vel_distr")
    end

    #
    # Next, create the actual particle instances based on the solver (full orbit
    # or GCA).
    #
    # Full orbit
    if params.solver == "full-orbit"
        particles = ParticleSoA(pos,
                                vel,
                                params.specie,
                                params.nsteps)
    elseif params.solver == "full-orbit-isf"
        particles = ParticleSoA(pos,
                                vel,
                                params.specie,
                                params.nsteps)
    elseif params.solver == "full-orbit-relativistic"
        particles = ParticleSoA(pos,
                                vel,
                                params.specie,
                                params.nsteps)
    # GCA
    elseif params.solver == "GCA"
        vparal = zeros(params.wp_part, numdims, params.npart)
        magneticMoment = zeros(params.wp_part, params.npart)
        for i = 1:params.npart
            (B⃗, E⃗), _ = gridinterp(mesh,
                                   params.interp,
                                   pos[:,i]
                                   )
            B = norm(B⃗)
            b̂ = B⃗/B
            v = norm(vel[:,i])
            vparal[:,i] = [vel[:,i] ⋅ b̂]
            vperp = √(v^2 - vparal[:,i]^2)
            mass = specieTable[params.specie[i], 1]
            magneticMoment[i] = [mass*vperp^2/(2B)]
        end
        particles = GCAParticleSoA(pos,
                                   vparal,
                                   magneticMoment,
                                   params.specie,
                                   params.nsteps
                                   )
    else
        error("Invalid solver-parameter: $params.solver")
    end
        
    #---------------------------------------------------------------------------
    # Construct Patch
    #---------------------------------------------------------------------------
    patch = Patch(mesh,
                  particles,
                  methodmap[params.solver],
                  methodmap[params.scheme],
                  methodmap[params.interp],
                  params.dt,
                  params.nsteps,
                  params.npart,
                  params.pbc
                  )

    #---------------------------------------------------------------------------
    # Construct Experiment
    #---------------------------------------------------------------------------
    exp = Experiment(params,
                     patch
                     )
    return exp
end # function tp_init
#|

function tp_init(
    params::Parameters,
    mesh_filename::String,
    bg_filename::String,
    tp_filename::String,
    )
    #
    numdims = 3
    #
    println("tp.jl: Initialising experiment...")

    #
    # Open tp-file
    #
    f = open(tp_filename)
    pos = zeros(params.wp_part, numdims, params.npart, params.nsteps)
    vel = zeros(params.wp_part, numdims, params.npart, params.nsteps)
    pos[1,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    pos[2,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    pos[3,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    vel[1,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    vel[2,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    vel[3,:,:] = mmap(f, Matrix{params.wp_part}, (params.npart, params.nsteps))
    close(f)

    #
    # Open mesh-file
    #
    f = open(mesh_filename)
    x = mmap(f, Vector{params.wp_snap}, params.nx)
    y = mmap(f, Vector{params.wp_snap}, params.ny)
    z = mmap(f, Vector{params.wp_snap}, params.nz)
    close(f)
    
    #
    # Open bg-file
    #
    meshsize = (params.nx, params.ny, params.nz)
    bField = zeros(params.wp_snap, 3, meshsize...)
    eField = zeros(params.wp_snap, 3, meshsize...)
    ∇B = zeros(params.wp_snap, 3, meshsize...)
    ∇b̂ = zeros(params.wp_snap, 3, 3, meshsize...)
    ∇ExB = zeros(params.wp_snap, 3, 3, meshsize...)
    f = open(bg_filename)
    for i = 1:3
        bField[i,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
    end
    for i = 1:3
        eField[i,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
    end
    for i = 1:3
        ∇B[i,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
    end
    for i = 1:3
        for j = 1:3
            ∇b̂[i,j,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
        end
    end
    for i = 1:3
        for j = 1:3
            ∇ExB[i,j,:,:,:] = mmap(f, Array{params.wp_snap, 3}, meshsize)
        end
    end
    close(f)
    
    #
    # Create Mesh, ParticleSoA, Patch and Experiment
    #
    mesh = Mesh(bField, eField, ∇B, ∇b̂, ∇ExB, x, y, z)
    particles = ParticleSoA(pos, vel, params.specie)
    patch = Patch(mesh,
                  particles,
                  methodmap[params.solver],
                  methodmap[params.scheme],
                  methodmap[params.interp],
                  params.dt,
                  params.nsteps,
                  params.npart,
                  params.pbc
                  )
    exp = Experiment(params, patch);
    #
    return exp, params
end # function tp_init
#|

function tp_softinit!(
    exp::Experiment
    )
end


function tp_run!(
    exp::Experiment
    )
    statement = string("tp.jl: Running simulation:\n",
                       "\tnpart:  $(exp.params.npart)\n",
                       "\tnsteps: $(exp.params.nsteps)\n",
                       "\tdt:     $(exp.params.dt)\n",
                       "\tNumber of iterations: ",
                       "$(exp.params.npart*exp.params.nsteps)"
                       )
    println(statement)
    run!(exp.patch)
end # function tp_run
#|

function tp_run(
    params::Parameters
    )
    exp = tp_init!(params)
    tp_run!(exp)
    return exp
end # function tp_run


function tp_save(
    exp    ::Experiment
    ;
    expname=nothing,
    expdir=nothing,
    )
    #
    # Check whether path for saving files is defined.
    #
    if expname == nothing
        if !isdefined(exp.params, :tp_expname)
            error(string("Experiment name `tp_expname` needed for saving.",
                         " Specify it in `Parameters` or give as argument to",
                         " `tp_saveExp`."))
        else
            expname = exp.params.tp_expname
        end
    end
    if expdir == nothing
        if !isdefined(exp.params, :tp_expdir)
            error(string("Experiment directory `tp_expdir` needed for saving.",
                         " Specify it in `Parameters` or give as argument to",
                         " `tp_saveExp`."))
        else
            expdir = exp.params.tp_expdir
        end
    end

    # To avoid magic numbers
    numdims = 3

    println("tp.jl: Saving experiment...")
    #
    # Construct filenames
    #
    basename = string(expdir, "/", expname)
    tp_filename = string(basename, ".tp")
    mesh_filename = string(basename, ".mesh")
    bg_filename = string(basename, ".bg")
    params_filename = string(basename, "_params", ".jl")
    # SUGGESTION: Add some parameter-file? Not sure yet what to include and how
    # to make it. It would be nice with a params-file or script with same name
    # that could be run to produce the corresponding .tp file.

    #
    # Write particle positions, velocity and life-status
    #
    # IMPROVEMENT: This method takes a up a lot of memory. Defining a
    # getArrayForSaving() function in the Particles module could be better.
    f = open(tp_filename, "w+")
    write(f, exp.patch.tp.pos[1,:,:])
    write(f, exp.patch.tp.pos[2,:,:])
    write(f, exp.patch.tp.pos[3,:,:])
    write(f, exp.patch.tp.vel[1,:,:])
    write(f, exp.patch.tp.vel[2,:,:])
    write(f, exp.patch.tp.vel[3,:,:])
    close(f)
    println("tp.jl: Wrote $tp_filename")

    #
    # Write mesh (NB! here I seperate the fields from the mesh)
    #
    f = open(mesh_filename, "w+")
    write(f, [exp.patch.mesh.xCoords[:];
              exp.patch.mesh.yCoords[:];
              exp.patch.mesh.zCoords[:];
              ]
          )
    close(f)
    println("tp.jl: Wrote $mesh_filename")

    #
    # Write background fields
    #
    f = open(bg_filename, "w+")
    write(f, exp.patch.mesh.bField[1,:,:,:])
    write(f, exp.patch.mesh.bField[2,:,:,:])
    write(f, exp.patch.mesh.bField[3,:,:,:])
    write(f, exp.patch.mesh.eField[1,:,:,:])
    write(f, exp.patch.mesh.eField[2,:,:,:])
    write(f, exp.patch.mesh.eField[3,:,:,:])
    write(f, exp.patch.mesh.∇B[1,:,:,:])
    write(f, exp.patch.mesh.∇B[2,:,:,:])
    write(f, exp.patch.mesh.∇B[3,:,:,:])
    for i = 1:3
        for j = 1:3
            write(f, exp.patch.mesh.∇b̂[i,j,:,:,:])
        end
    end
    for i = 1:3
        for j = 1:3
            write(f, exp.patch.mesh.∇ExB[i,j,:,:,:])
        end
    end
    close(f)
    println("tp.jl: Wrote $bg_filename")

    #
    # Write parameters
    #
    paramsstring = "using TraceParticles\nparams = Parameters(\n"
    for p in fieldnames(Parameters)
        if isdefined(exp.params, p)
            if typeof(getfield(exp.params, p)) == String
                value = "\"$(getfield(exp.params, p))\""
            else
                value = "$(getfield(exp.params, p))"
            end
            spaces = " "^(11 - length("$p"))
            paramsstring = string(paramsstring,
                                  "\t$p", spaces, "= $value,\n")
        end
    end
    paramsstring = string(paramsstring, ")")
    f = open(params_filename, "w+")
    write(f, paramsstring)
    close(f)
    println("tp.jl: Wrote $(params_filename)")

    #
end # function tp_saveExp


function tp_load(
    params ::Parameters,
    ;
    expdir ::String=params.tp_expdir,
    expname::String=params.tp_expname,
    )
    
    basename = string(expdir, "/", expname)
    tp_filename = string(basename, ".tp")
    mesh_filename = string(basename, ".mesh")
    bg_filename = string(basename, ".bg")
    #params_filename = string(basename, ".jl")

    #
    # Get parameters
    #
    #f = open(params_filename)
    #contents = read(f, String)
    #commands = Meta.parse("begin $contents end").args
    #println("upd")
    #for line in commands
    #    eval(line)
    #end

    #
    # Create experiment and 
    #
    exp = tp_initfromfile(params, mesh_filename, bg_filename, tp_filename);
    #
    return exp
end # function tp_loadExp


function tp_set_dt!(
    exp::Experiment,
    dt ::Real,
    )
    exp.patch.dt = dt
end


#-------------------------------------------------------------------------------
end # module tp


