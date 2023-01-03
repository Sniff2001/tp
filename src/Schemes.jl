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
    nextPos = @. pos + vel * dt
    nextVel = @. vel + acc * dt
    return nextPos, nextVel
end

function eulerCromer(pos, vel, acc, dt)
    nextVel = @. vel + acc * dt
    nextPos = @. pos + nextVel * dt
    return nextPos, nextVel
end

end
