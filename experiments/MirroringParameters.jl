module MirroringParameters

using Utilities
using Particles
using WorkingPrecision
using TraceParticles
using Test
using Interpolations_tp: trilinear_ip, quadratic_bspline, cubic_bspline
using Meshes:           Mesh
using Solvers          
using Schemes
using Interpolations_tp
using Utilities
using LinearAlgebra

export mirroringParameters

function mirroringParameters(n, dt; interpolator="trilinear")
    numparticles = 1  # Number of particles to simulate
    numdims = 3
    species = 4*ones(wpInt, numparticles) # Specifies the species of the particles
    mass = specieTable[species[1], 1]
    charge = specieTable[species[1], 2]
    #...............................................
    # MAGNETIC FIELD PARAMETERS
    B0 = 30.0 # Magnetic field strenth parameter
    L = 0.4 # Mirroring length
    vel0 = [0.0, 0.1, 0.1]
    rL = mass*√(vel0[1]^2 + vel0[2]^2)/(charge*B0)
    # Lower bounds of the three spatial axes
    a = 2.1L
    xi0 = (-a, -a, -a)
    # Upper bound of the three spatial axes
    xif = (a, a, a)
    xx, yy, zz, dx, dy, dz = Utilities.createaxes(xi0, xif, n)
    Ex = 0.0
    Ey = 0.0
    Bfield = zeros(wpFloat, numdims, n[1], n[2], n[3])
    Efield = zeros(size(Bfield))
    Efield[1,:,:,:] .= Ex
    Efield[2,:,:,:] .= Ey
    #Efield[1,:,1:30,:] .= 0.0
    Utilities.discretise!(Bfield, xx, yy, zz, Utilities.mirroringfield, B0, L)
    Bx = Bfield[1,:,:,:]
    By = Bfield[2,:,:,:]
    Bz = Bfield[3,:,:,:]
    Ex = Efield[1,:,:,:]
    Ey = Efield[2,:,:,:]
    Ez = Efield[3,:,:,:]
    time = 30
    # Set simulation parameters
    params = Parameters(
    dt=dt,
    nsteps=round(Int, time/dt),
    npart=1,
    tp_expname="mirroring$(n[1])U$(dt)",
    tp_expdir="C:/Users/ixyva/data",
    solver="full-orbit",
    scheme="RK4",
    interp=interpolator,
    wp_part=Float32,
    wp_snap=Float32,
    pos_distr="point",
    vel_distr="point",
    posxbounds=[-rL],
    posybounds=[0.0],
    poszbounds=[0.0],
    velxbounds=[vel0[1]],
    velybounds=[vel0[2]],
    velzbounds=[vel0[3]],
    x=xx,
    y=yy,
    z=zz,
    bx=Bx,
    by=By,
    bz=Bz,
    ex=Ex,
    ey=Ey,
    ez=Ez,
    nx=n[1],
    ny=n[2],
    nz=n[3],
    bg_input="manual",
    specie=4*ones(Int64, 1),
    pbc=(true, true, true)
    )
    return params
end

if abspath(PROGRAM_FILE) == @__FILE__
    using TPplots
    import PyPlot
    const plt = PyPlot
    plt.pygui(true)

    n = (100, 100, 100)
    dt = 1.e-3        # Time step [s]
    time = 30.0 #n=100   # End time of simulation [s]
    #time = 2.3 #n=10   # End time of simulation [s] 
    params = mirroringParameters(n, dt)
    # Initialise experiment
    exp = tp_init!(params)
    # Run experiment
    tp_run!(exp)
    #plot(exp.patch, "x")
    plot3D(exp.patch)
    plt.show()

    B0 = 30.0 # Magnetic field strenth parameter
    numparticles = 1
    species = 4*ones(wpInt, numparticles) # Specifies the species of the particles
    L = 0.4
    a = 2.1L
    xi0 = (-a, -a, -a)
    # Upper bound of the three spatial axes
    xif = (a, a, a)
    xx, yy, zz, dx, dy, dz = Utilities.createaxes(xi0, xif, n)
    numdims = 3
    Ex = 0.0
    Ey = 0.0
    Bfield = zeros(wpFloat, numdims, n[1], n[2], n[3])
    Efield = zeros(size(Bfield))
    Efield[1,:,:,:] .= Ex
    Efield[2,:,:,:] .= Ey
    #Efield[1,:,1:30,:] .= 0.0
    Utilities.discretise!(Bfield, xx, yy, zz, Utilities.mirroringfield, B0, L)
    mass = specieTable[species[1], 1]
    charge = specieTable[species[1], 2]
    vel0 = [0.0, 0.1, 0.1]
    rL = mass*√(vel0[1]^2 + vel0[2]^2)/(charge*B0)
    pbc = (true, true, true)

    #interp = trilinear_ip(Bfield, Efield, xx, yy, zz, pbc)
    
    #B⃗ = interp[1](-rL, 0, 0)
    #E⃗ = interp[2](-rL, 0, 0)

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

    function z(t, μ, mass, B0, L, A, ϕ)
        ω = √(2μ*B0/(mass*L^2))
        return @. A*sin(ω*t + ϕ)
    end
    ϕ = 0
    A = zmax
    times = collect(range(0.0, step=dt, length=round(Int, time/dt)+1))
    z_anal = z(times, μ[1], mass, B0, L, A, ϕ)
    z_FO = exp.patch.tp.pos[3,1,:]
    rmse_FO = √(sum((z_anal .- z_FO).^2)/round(Int, time/dt))
    @testset verbose = true "Full orbit: RK4" begin
        @test isapprox(rmse_FO, 0.0, atol=0.001)
    end # testset GCA: Euler
end

end