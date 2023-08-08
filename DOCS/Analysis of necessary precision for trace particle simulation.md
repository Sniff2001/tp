## Method

The script takes a high precision value as the "true" value and compares it to each simulation for different precisions. For the first two comparisons we will check grid size and timestep. Here it is assumed that the perfect true solution is when the grid size and timestep approaches infinity. However, as this isn't computably feasible we decide to use the true value as a relatively small value.

### Interpolations
Here is a brief overview of how to use Interpolations.jl in the form of how Interpolations_tp.jl is made:

```jl
using Interpolations
itp = interpolate((xx,yy,zz), tensorfield, Gridded(Linear()))
val = itp(x,y,z)
```
This performs a trilinear interpolation, however, it seems to interpolate the entire field, so it is slower than Interpolations_tp.jl and therefore does the entire interpolation ~4x slower. If you want to interpolate for a smaller dataset you probably have to choose only the field that's relevant.
The second line; `interpolate()` takes in the grid size as a tuple of arrays (optional, although seems to only be for B-splines where you scale it afterwards), the field to interpolate and the interpolation methods.
The third line evaluates the interpolation at a given position. The following interpolations are available `Gridded()` and `BSpline()`. What you include within them determines their dimension:
- `Gridded(Constant())` - Nearest neighbour
- `Gridded(Linear())` - Trilinear interpolation for 3D grid, bilinear interpolation for 2D grid
- `BSpline(Constant())` - B-spline of degree 0 (equivalent to nearest neighbour)
- `BSpline(Linear())` - B-spline of degree 1 (equivalent to trilinear interpolation)
- `BSpline(Quadratic())` - B-spline of degree 2
- `BSpline(Cubic())` - B-spline of degree 3

Since B-splines of degree 2 or higher uses more points than the closest, boundary conditions are needed. Luckily, Interpolations.jl got you covered. `Flat()`, `Line()`, `Free()`, `Periodic()` and `Reflect()` are the options for boundary conditions. It is also required to choose where the boundary condition is applied, normally on the boundary point you'd use `OnGrid()`, but if the point represents the middle of the grid then `OnCell()`. The boundary conditions are used inside the interpolation: i.e. `BSpline(Cubic(Periodic(OnGrid())))`.

Summary from here and trial and error:
https://github.com/JuliaMath/Interpolations.jl/blob/master/docs/src/control.md

## Results & Discussion

### Rewrite
The rewrite of the PlasmoidParameters.jl and DipoleParameters.jl seems to give the correct plot, the base plot being:
![[Pasted image 20230628141228.png|400x300]]
![[Pasted image 20230628140622.png|400x300]]

These will be relevant for comparison with future plots.

### The impact of grid size
![[Pasted image 20230628155730.png]]
(Using a loglog plot here for better visibility)
Here the grid sizes vary from 2 to 1e3, where 1e3 is chosen as the value the errors are compared against. Due to it being an approximate linear plot in loglog, it might be fair to assume that the general formula is an exponential $\epsilon(x) = x^a + 10^b$ which reduces to $\log(\epsilon(x)) = a\log(x) + b$ where $\log(x)$ is the logarithm of 10. The grid size seems to have a minor impact on the results, showing errors within a 1% margin at the standard of 100x100x2. However, decreasing the grid size below 100 doesn't seem to heavily affect it too badly. Even at 10x10x2 (plot below) it shows relatively the same path. Only notable difference is particle 2 and 3.
![[Pasted image 20230629105055.png]]


### The impact of timestep
![[Pasted image 20230628143902.png|700x450]]
The variance in solution from timestep is plotted above. The x-axis showing the different timestep sizes while the relative error on the y-axis in percentage. The dataset is 60 logarithmic values (removed duplicates, might be smaller) from 1e-5 to 1. Here the timestep 1e-5 is chosen as the value the error is compared against.
It seems that 1e-3 is a relatively sufficient choice of timestep, although particle 1 reaches errors up to 10% of the 1e-5 path. A better choice might be 5e-4, which is the closest point before the spikes start appearing. The values lower than this threshold only have a relative error of around 0.5%. It is also interesting to note when the jumps on the graph happens. One might notice that particle 1 needs a lower timestep than the others to be accurate. Particle 4 barely needs a timestep to be accurate and is sufficient for the 1e-3. The immense curl and stronger magnetic field particle 1, 2 and 3 encounters might be the reason for this and will have a higher effect on the final position of the particle. It is also interesting that the graph seems to show that it is better to choose a timestep of 0.1 over 0.01. Saving the most interesting for last is why the graph totally flattens out at high values. Here are a few plots that show the paths at chosen timesteps of: 1e-5, 1e-2, 5e-2 and 1.

![[Pasted image 20230628153642.png]]
plot of 1e-5, the most accurate plot.

![[Pasted image 20230628153803.png]]
plot of 1e-2

![[Pasted image 20230628153913.png]]
plot of 5e-2

![[Pasted image 20230628154028.png]]
plot of 1

Interestingly enough, it seems as if everything above 1e-2 doesn't have enough time to react compared to 1e-5. This is because the particles are immediately killed from step 1. These solutions might therefore just coincidentally be a closer approximation than the ones that actually have time. It also makes sense for them to flatten out as none higher than a certain value will have the time to react at all. It is still easy to see from the plot of 1e-2 seems to be able to plot particle 4 with minor deviations from 1e-5. However, the other particles seem to struggle the closer they get to the stronger magnetic field.

### Static vs interdynamic field
For a visual clue, let us visualise how a static field impacts the path for default values (n = 100x100x2 and dt = 0.001):
![[Pasted image 20230629110134.png]]

It is quite clear that the static field doesn't work for the current initial values. How does this compare to timestep and grid size? We can find the same plots for a static field, let us first visualise what the plot even looks like at the actual value that we compare the error against:
![[Pasted image 20230629110825.png]]
changed grid size (remind me to ask eilif if there's a fix for the quiver, and if not, implement it)

![[Pasted image 20230629111029.png]]
changed timestep

It is very clear that this is a timestep issue from eyeing it, which makes sense as timestep allows it to check the magnetic field more often. Since this practically works like interdynamic field by checking the magnetic field more often. However, which takes more time one should consider which method is better.

![[Pasted image 20230629111844.png]]

![[Pasted image 20230629113023.png]]

In the dt I changed it to loglog as a very rogue particle 4 suddenly decided to jump to >100% error.

It is easy to see that at dt = 1e-3 the error is quite significant at around 10%. The grid size also seems to give higher errors for more grid sizes than with interdynamic. What if we choose the actual value to be the interdynamic one?
![[Pasted image 20230629114914.png]]

![[Pasted image 20230629114741.png]]

But it seems pretty clear from this example that the interdynamic is simply better, static field seems to be worse. Again, could depend on if changing grid size and timestep saves more time than changing to interdynamic, however, these must be rare cases.
The grid size on the other hand seems to converge at around 100x100x2, seemingly unable to be beneficial for higher values. Showing that the error might always have a fixed value from the interdynamic for smaller grid sizes. It is unclear why that might be the case.

### Interpolations

1st degree B-spline interpolation is slower than trilinear interpolation directly by a few tests. Despite the Interpolations.jl being slower finding an interpolation, it might be beneficial to use if the magnetic/electric field is unchanging, static over time, since the actual evaluation is much much faster than calling Interpolations_tp.jl for each step.

Benchmarks for the plasmoid experiment at default parameters:
![[image.png]]
^ Using trilinear interpolation

![[image 1.png]]
^ Using trilinear from Interpolations.jl

Benchmarks for dipole experiment at default parameters:
![[Pasted image 20230727135624.png]]
^ Using trilinear interpolation

![[Pasted image 20230727135851.png]]
^ Using trilinear from Interpolations.jl

Not only is it faster, but it also uses less memory. The two cases I can imagine the old interpolations library is more beneficial is as said with a changing field and for small simulations (in terms of steps).

Moving forward, we will check if you can neglect some resolution in favour of a higher order interpolation.

![[Pasted image 20230727155926.png]]
Grid size of 20x20x2 (from left: trilinear, quadratic bspline, cubic bspline)

Above is a quick image of what we are looking at. Below we will be comparing the relative error compared to a highly resolved grid.

![[Pasted image 20230802145820.png]]
The different interpolators relative error compared to a trilinear n=1000x1000x2 and dt=1e-5

As seen above the interpolators tend to converge toward a fixed error. This is possibly the limit for how much you can improve the error without changing the time step. Although, it can be seen that the B-splines converge much smoother than trilinear interpolation. There is therefore less error between different grid sizes using B-splines when approaching convergence contrary to trilinear which is oscillating heavily around convergence using the same grid size. The explanation could be that B-splines keep a relatively similar curve and better at preserving these compared to trilinear which receives a greater curve the more points there are.

![[Pasted image 20230804092243.png]]

Here we compare the quadratic and cubic B-spline against trilinear of the same time step. Strangely enough, despite being similar to the time step at 1e-5, particle 4 shows the B-splines reach a more value at a lower grid size. As for the rest, the trilinear oscillates heavily as before.

What is interesting about this is that particle 4 is exposed to the least gradient of the magnetic field. This could indicate that B-splines are better at predicting small changes in the magnetic field. 

Let us compare the magnetic field calculated against an analytical solution to see which has the least error.
![[image 2.png]]
This shows the plasmoid at default parameters (but I only let it run for 10 steps) (also it technically starts at 1 but ignore that) The error is found by taking the norm of the magnetic field vector and comparing these two values.

This shows a much clearer error in trilinear compared to the splines, while not much, it is something. It is most notable for particle 3 where the relative error almost doubles. To be frank, running this code takes forever, possibly because it is not optimised. So I let one run go to 500 steps max:
![[Pasted image 20230808142340.png]]
While the errors are less than a percentage, it is worth noting that these accumulate over time. Again, trilinear is shown to be very unstable unlike the B-splines that keep a smoother error over time.

Lastly I wanted to check if maybe the energy conservation is better with the interpolators:
![[Pasted image 20230804140732.png]]
Trilinear (100000 particles)

![[Pasted image 20230804112318.png]]
Cubic B-spline (100000 particles)

Although, while different, it is practically the same. Meaning the total energy on a random set of particles at different initial conditions barely change overall. This is expected considering the interpolation errors are basically 0.1% from the analytical solution.

## Appendix
The requested loglog plot of timestep:

![[Pasted image 20230629120747.png]]

![[Pasted image 20230731114921.png]]
![[Pasted image 20230731113907.png]]
Quadratic bsplines compared to trilinear standard at n=(1000, 1000, 2), dt=1e-5

![[Pasted image 20230731115810.png]]
![[Pasted image 20230731113600.png]]
Cubic bsplines compared to trilinear standard (same as above)





