module LineBroadenings

import SpecialFunctions: erfcx
import PhysicalConstants.CODATA2018: k_B as 𝑘, c_0 as 𝑐
using Unitful
import Unitful: Mass, AbsoluteScaleTemperature, Pressure

export Lorentz, Gaussian, Voigt, pdf, peak, fwhm
include("profiles.jl")

export profile_doppler, profile_pressure, profile_voigt, fwhm_doppler, fwhm_pressure
include("spectra.jl")

end
