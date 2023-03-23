abstract type Profile end

Broadcast.broadcastable(d::Profile) = Ref(d)

"""
    pdf(d::Profile, x::Real)

probability density function
"""
pdf

include("gaussian.jl")
include("lorentz.jl")
include("voigt.jl")
