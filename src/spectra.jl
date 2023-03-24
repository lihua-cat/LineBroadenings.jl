"""
    profile_pressure(ν; ν0, νp)

line profile of pressure broadening (Lorentz profile).

# Arguments
- `ν`: target wavelength or wavenumber
- `ν0`: central wavelength or wavenumber
- `νp`: FWHM of pressure broadening
"""
function profile_pressure(ν; ν0, νp)
    u = unit(ν)
    ν0, νp = uconvert.(u, (ν0, νp))
    ν, ν0, νp = ustrip.((ν, ν0, νp))
    d_l = Lorentz(ν0, νp / W_L)
    pdf(d_l, ν) / u
end

"""
    profile_doppler(ν; ν0, νd)

line profile of Doppler broadening (Gaussian profile).

# Arguments
- `ν`: target wavelength or wavenumber
- `ν0`: central wavelength or wavenumber
- `νd`: FWHM of Doppler broadening
"""
function profile_doppler(ν; ν0, νd)
    u = unit(ν)
    ν0, νd = uconvert.(u, (ν0, νd))
    ν, ν0, νd = ustrip.((ν, ν0, νd))
    d_g = Gaussian(ν0, νd / W_G)
    pdf(d_g, ν) / u
end

"""
    profile_voigt(ν; ν0, νd, νp)

line profile of combined Doppler and pressure broadenings (Voigt profile).

# Arguments
- `ν`: target wavelength or wavenumber
- `ν0`: central wavelength or wavenumber
- `νd`: FWHM of Doppler broadening
- `νp`: FWHM of pressure broadening
"""
function profile_voigt(ν; ν0, νd, νp)
    u = unit(ν)
    ν0, νd, νp = uconvert.(u, (ν0, νd, νp))
    ν, ν0, νp, νd = ustrip.((ν, ν0, νp, νd))
    d_v = Voigt(ν0, νd / W_G, νp / W_L)
    pdf(d_v, ν) / u
end

@doc raw"""
    fwhm_pressure(P, γ)

FWHM of pressure broadening.

```math
\nu_p = 2 \sum P \cdot γ
```

# Arguments
- `P::Pressure`: pressure
- `γ`: pressure broadening coefficients of HWHM
"""
fwhm_pressure(P::Pressure, γ) = 2 * sum(P .* γ)

@doc raw"""

    fwhm_doppler(ν0, M, T)

FWHM of Doppler broadening.

```math
\nu_d = \sqrt{\frac{8kT\ln{2}}{Mc^2}}\nu_0
```

# Arguments
- `ν0`: central wavelength or wavenumber
- `M`: atomic weight
- `T`: temperature in Kelvin
"""
function fwhm_doppler(ν0, M::Mass, T::AbsoluteScaleTemperature)
    uconvert(unit(ν0), √(8𝑘 * T * log(2) / M / 𝑐^2) * ν0)
end

"approx fwhm of Voigt profile"
fwhm_voigt_approx(νd, νc) = uconvert(unit(νd), 0.5346νc + √(0.2166νc^2 + νd^2))
