module LineBroadenings

import SpecialFunctions: erfcx

export Lorentz, Gaussian, Voigt, pdf, peak, fwhm
include("profiles.jl")

end
