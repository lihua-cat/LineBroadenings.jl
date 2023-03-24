abstract type Profile end

Broadcast.broadcastable(d::Profile) = Ref(d)

"""
    pdf(d::Profile, x::Real)

probability density function.
"""
pdf

"""
    peak(d::Profile)

maximum of the pdf of `d`.
"""
peak

"""
    fwhm(d::Profile)

The full width at half maximum (FWHM) of `d`.
"""
fwhm

include("gaussian.jl")
include("lorentz.jl")
include("voigt.jl")
