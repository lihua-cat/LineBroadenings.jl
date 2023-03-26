##
using GLMakie
using Unitful
import PhysicalConstants.CODATA2018: c_0 as 𝑐
using LineBroadenings
##
μ, γ, σ = 0, 20 * √(2 * log(2)), 20
fl = Lorentz(μ, γ)
fg = Gaussian(μ, σ)
fv = Voigt(μ, σ, γ)
x = LinRange(-250, 250, 100000)
##
let
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, x, pdf.(fl, x), label = "Lorentz")
    lines!(ax, x, pdf.(fg, x), label = "Gaussian")
    lines!(ax, x, pdf.(fv, x), label = "Voigt")
    axislegend()
    fig
    # save("fig1.png", fig, px_per_unit = 2)
end
##
M = 126.9u"u"
k34 = 7603.14u"cm^-1"

kx = collect(7603.05:0.00001:7603.23)u"cm^-1"
T, P = 150 * u"K", 10 * u"Torr"

fwhm_d = fwhm_doppler(k34, M, T)
fwhm_p = fwhm_pressure(P, 5.0u"MHz/Torr") / 𝑐 |> u"cm^-1"
##
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Wavenumber (cm⁻¹)", ylabel="intensity distribution(a. u.)", title="line broadenings of I¹²⁷ (T=$T, P=$P)")
    offset = 7603
    kxx = ustrip.(u"cm^-1", kx) .- offset
    lines!(ax, kxx, ustrip.(profile_pressure.(kx; ν0=k34, νp=fwhm_p)), label="Pressure")
    lines!(ax, kxx, ustrip.(profile_doppler.(kx; ν0=k34, νd=fwhm_d)), label="Doppler")
    lines!(ax, kxx, ustrip.(profile_voigt.(kx; ν0=k34, νd=fwhm_d, νp=fwhm_p)), label="Voigt")
    ax.xtickformat = xs -> ["$(x + offset)" for x in xs]
    axislegend()
    fig
    # save("fig2.png", fig, px_per_unit = 2)
end
