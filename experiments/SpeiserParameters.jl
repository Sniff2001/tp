module SpeiserParameters

using Utilities
using Particles
using WorkingPrecision
using TraceParticles

export speiserParameters

η = 0.025  # guide field?
d = 1e-4 # Current sheet width
b = 1e-2  # Characteristic field strength

function speiserParameters(n, dt; interpolator="trilinear")
    mass = specieTable[species[1], 1]
    charge = specieTable[species[1], 2]
    t_eject = π*mass/(charge*η*b)
    time = 1.1*t_eject
    v0 = 1e7 
    a = 1e-2 * v0*b
    Ez = -a
    qMm = 20.0
    M = qMm*mass/charge
    rL = mass/(charge*M)
    a = 1.5
    L0 = 1e4 # m
    xi0 = (-3e-6*L0, -0.2*L0, -0.2*L0)
    xif = (3e-6*L0, 0.2*L0, 0.2*L0)
    xx, yy, zz, dx, dy, dz = Utilities.createaxes(xi0, xif, n)
    Bfield = zeros(Float64, 3, n[1], n[2], n[3])
    Efield = zeros(size(Bfield))
    Efield[3, :,:,:] .= Ez
    Utilities.discretise!(Bfield, xx, yy, zz, speiserBfield)
    Bx = Bfield[1,:,:,:]
    By = Bfield[2,:,:,:]
    Bz = Bfield[3,:,:,:]
    Ex = Efield[1,:,:,:]
    Ey = Efield[2,:,:,:]
    Ez = Efield[3,:,:,:]
    # Set simulation parameters
    params = Parameters(
    dt=dt,
    nsteps=round(Int, time/dt),
    npart=1,
    tp_expname="dipole$(n[1])U$(dt)",
    tp_expdir="C:/Users/ixyva/data",
    solver="full-orbit",
    scheme="RK4",
    interp=interpolator,
    wp_part=Float32,
    wp_snap=Float32,
    pos_distr="point",
    vel_distr="point",
    posxbounds=[1e-6*L0],
    posybounds=[1e-10*L0],
    poszbounds=[1e-10*L0],
    velxbounds=[0.0],
    velybounds=[0.0],
    velzbounds=[0.0],
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

function speiserBfield(
    x::wpFloat,
    y::wpFloat,
    z::wpFloat
    )
    return [η, -x/d, 0.0]*b
end

if abspath(PROGRAM_FILE) == @__FILE__
    using TPplots
    import PyPlot
    const plt = PyPlot
    plt.pygui(true)

    n = (100, 100, 100)
    dt = 0.2e-8
    numparticles = 1
    species = 4*ones(wpInt, numparticles)
    time = 100
    params = speiserParameters(n, dt)
    # Initialise experiment
    exp = tp_init!(params)
    # Run experiment
    tp_run!(exp)
    #plot(exp.patch, "y")
    plot3D(exp.patch)
    plt.show()
end

end