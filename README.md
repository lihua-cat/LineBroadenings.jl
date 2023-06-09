# LineBroadenings

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lihua-cat.github.io/LineBroadenings.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lihua-cat.github.io/LineBroadenings.jl/dev/)
[![Build Status](https://github.com/lihua-cat/LineBroadenings.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lihua-cat/LineBroadenings.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/lihua-cat/LineBroadenings.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/lihua-cat/LineBroadenings.jl)


## Example

```julia
using GLMakie
using LineBroadenings

μ, γ, σ = 0, 20 * √(2 * log(2)), 20
fl = Lorentz(μ, γ)
fg = Gaussian(μ, σ)
fv = Voigt(μ, σ, γ)
x = LinRange(-250, 250, 100000)

let
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, x, pdf.(fl, x), label = "Lorentz")
    lines!(ax, x, pdf.(fg, x), label = "Gaussian")
    lines!(ax, x, pdf.(fv, x), label = "Voigt")
    axislegend()
    fig
end
```

![fig1](/example/fig1.png)

```julia
M = 126.9u"u"
k34 = 7603.14u"cm^-1"

kx = collect(7603.05:0.00001:7603.23)u"cm^-1"
T, P = 150 * u"K", 10 * u"Torr"

fwhm_d = fwhm_doppler(k34, M, T)
fwhm_p = fwhm_pressure(P, 5.0u"MHz/Torr") / 𝑐 |> u"cm^-1"
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
end
```

![fig2](example/fig2.png)