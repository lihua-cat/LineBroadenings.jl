import PhysicalConstants.CODATA2018: c_0 as 𝑐

using LineBroadenings
using Unitful
using Test

@testset verbose=true "line profile functions" begin
    μ, γ, σ = rand(3)
    μ = (μ - 0.5) * 10
    fl = Lorentz(μ, γ)
    fg = Gaussian(μ, σ)
    fv = Voigt(μ, σ, γ)
    x = (rand(10_000_000) .- 0.5) * 1e3
    sort!(x)
    pl = pdf.(fl, x)
    pg = pdf.(fg, x)
    pv = pdf.(fv, x)
    @testset "constructor" begin
        @test_throws DomainError Lorentz(μ, -γ)
        @test_throws DomainError Gaussian(μ, -σ)
        @test_throws DomainError Voigt(μ, -σ, γ)
        @test_throws DomainError Voigt(μ, σ, -γ)
        @test_throws DomainError Voigt(μ, -σ, -γ)
        @test Voigt(μ, 0, γ) == Lorentz(μ, γ)
        @test Voigt(μ, σ, 0) == Gaussian(μ, σ)
        @test Lorentz(1, 1.0) == Lorentz(1, 1) == Lorentz(1.0, 1.0)
        @test Gaussian(1.0, 1) == Gaussian(1, 1) == Gaussian(1.0, 1.0)
        @test Voigt(0.0, 1.0, 1) == Voigt(0, 1, 1) == Voigt(0.0, 1.0, 1.0)
    end
    @testset "zero width" begin
        @test isinf(pdf(Lorentz(μ, 0), μ))
        @test isinf(pdf(Gaussian(μ, 0), μ))
        @test isinf(pdf(Voigt(μ, 0, 0), μ))
        @test iszero(pdf(Lorentz(μ, 0), μ + rand()))
        @test iszero(pdf(Gaussian(μ, 0), μ + rand()))
        @test iszero(pdf(Voigt(μ, 0, 0), μ + rand()))
    end
    @testset "peak" begin
        @test peak(fl) == 1 / π / γ
        @test peak(fg) == 1 / σ / √(2π)
        @test peak(fv) < peak(fl) && peak(fv) < peak(fg)
        @test maximum(pl)≈peak(fl) rtol=1e-3
        @test maximum(pg)≈peak(fg) rtol=1e-3
        @test maximum(pv)≈peak(fv) rtol=1e-3
    end
    @testset "center" begin
        @test x[findmax(pl)[2]]≈μ rtol=5e-3
        @test x[findmax(pg)[2]]≈μ rtol=5e-3
        @test x[findmax(pv)[2]]≈μ rtol=5e-3
    end
    @testset "normalized" begin
        @test sum(pl[2:end] .* diff(x))≈1 atol=1e-2
        @test sum(pg[2:end] .* diff(x))≈1 atol=1e-2
        @test sum(pv[2:end] .* diff(x))≈1 atol=1e-2
    end
    @testset "fwhm" begin
        x1 = view(x, x .< μ)
        x2 = view(x, x .>= μ)
        pl1 = view(pl, x .< μ)
        pl2 = view(pl, x .>= μ)
        fwhm_l = x2[(findmin(abs.(pl2 .- peak(fl) / 2)))[2]] -
                 x1[(findmin(abs.(pl1 .- peak(fl) / 2)))[2]]
        pg1 = view(pg, x .< μ)
        pg2 = view(pg, x .>= μ)
        fwhm_g = x2[(findmin(abs.(pg2 .- peak(fg) / 2)))[2]] -
                 x1[(findmin(abs.(pg1 .- peak(fg) / 2)))[2]]
        pv1 = view(pv, x .< μ)
        pv2 = view(pv, x .>= μ)
        fwhm_v = x2[(findmin(abs.(pv2 .- peak(fv) / 2)))[2]] -
                 x1[(findmin(abs.(pv1 .- peak(fv) / 2)))[2]]

        @test fwhm_l≈fwhm(fl) rtol=1e-2
        @test fwhm_g≈fwhm(fg) rtol=1e-2
        @test fwhm_v≈LineBroadenings.fwhm_approx(fv) rtol=1e-2
        @test fwhm_v > fwhm(fl) && fwhm_v > fwhm(fg)
    end
end;

@testset verbose=true "line broadenings" begin
    M = 126.9u"u"
    k34 = 7603.14u"cm^-1"
    T, P = 150u"K", 10u"Torr"
    fwhm_d = fwhm_doppler(k34, M, T)
    fwhm_p = fwhm_pressure(P, 5.0u"MHz/Torr") / 𝑐 |> u"cm^-1"
    @testset "FWHM" begin
        @test fwhm_d ≈ 0.005920548179262646u"cm^-1"
        @test fwhm_p ≈ 0.0033356409519815205u"cm^-1"
    end
    @testset "profile" begin
        dk = 1e-7 * u"cm^-1"
        kx = (k34 - 1u"cm^-1"):dk:(k34 + 1u"cm^-1")
        pd = profile_doppler.(kx; ν0 = k34, νd = fwhm_d)
        pp = profile_pressure.(kx; ν0 = k34, νp = fwhm_p)
        pv = profile_voigt.(kx; ν0 = k34, νd = fwhm_d, νp = fwhm_p)
        σ = fwhm_d / (2 * √(2 * log(2)))
        γ = fwhm_p / 2
        @test maximum(pd) == 1 / σ / √(2π)
        @test maximum(pp) == 1 / π / γ
        @test maximum(pv) < maximum(pp) && maximum(pv) < maximum(pd)
        @test sum(pd .* dk)≈1 atol=0.01
        @test sum(pp .* dk)≈1 atol=0.01
        @test sum(pv .* dk)≈1 atol=0.01
    end
end;
