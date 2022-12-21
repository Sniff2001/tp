#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 Schemes.jl
#
#-------------------------------------------------------------------------------
# Module containing the various schemes for solving differential equations
#-------------------------------------------------------------------------------

module Schemes

function euler(pos, vel, acc, dt)
    @. pos = pos + vel * dt
    @. vel = vel + acc * dt
    return pos, vel
end

function eulerCromer(pos, vel, acc, dt)
    @. vel = vel + acc * dt
    @. pos = pos + vel * dt
    return pos, vel
end

end
