- [x] change dt
- [x] plasmoid.jl to datafile
- [x] define pos 3-vec
- [x] visualise differences in grid-size

benchmarktools -> @benchmark
@time func() -> prints time

-----

calculate interpolation in patches.jl -> run!()
add new method to eomLorentzforce() that doesn't require interpolator
???
profit.

- Union of function and ScaledInterpolation in patches.jl and particles.jl
- (`func4(x::Union{Function, Real}) = println(x)`)
- New method of fullorbit that uses the other method of eomLorentzforce
- (and add dictionary keys for new interpolations)