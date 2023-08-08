module DipoleParameters

using Utilities
using Particles
using WorkingPrecision
using TraceParticles

export dipoleParameters

function dipoleParameters(n, dt; interpolator="trilinear")
    numparticles = 1
    time = 100
    species = 4*ones(wpInt, numparticles)
    mass = specieTable[species[1], 1]
    charge = specieTable[species[1], 2]
    qMm = 20.0
    M = qMm*mass/charge
    rL = mass/(charge*M)
    a = 1.5
    xi0 = (-a, -a, -a)
    xif = (a, a, a)
    xx, yy, zz, dx, dy, dz = Utilities.createaxes(xi0, xif, n)
    Bfield = zeros(Float64, 3, n[1], n[2], n[3])
    Efield = zeros(size(Bfield))
    Utilities.discretise!(Bfield, xx, yy, zz, Utilities.dipolefield, M)
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
    posxbounds=[1-rL],
    posybounds=[0.0],
    poszbounds=[0.0],
    velxbounds=[0.0],
    velybounds=[1.0],
    velzbounds=[0.5],
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

    n = (256, 256, 256)
    dt = 5.e-4 
    params = dipoleParameters(n, dt, "trilinear")
    # Initialise experiment
    exp = tp_init!(params)
    # Run experiment
    tp_run!(exp)
    #plot(exp.patch, "y")
    plot3D(exp.patch)
    plt.show()
end

end