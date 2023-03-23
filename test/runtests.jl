using LineBroadenings
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
        @test peak(fl) == 1 / π / γ && isapprox(maximum(pl), peak(fl), rtol = 1e-3)
        @test peak(fg) == 1 / σ / √(2π) && isapprox(maximum(pg), peak(fg), rtol = 1e-3)
        @test isapprox(maximum(pv), peak(fv), rtol = 1e-3) &&
              peak(fv) < peak(fl) && peak(fv) < peak(fg)
    end
    @testset "center" begin
        @test isapprox(x[findmax(pl)[2]], μ, rtol = 5e-3)
        @test isapprox(x[findmax(pg)[2]], μ, rtol = 5e-3)
        @test isapprox(x[findmax(pv)[2]], μ, rtol = 5e-3)
    end
    @testset "normalized" begin
        @test isapprox(sum(pl[2:end] .* diff(x)), 1, rtol = 1e-2)
        @test isapprox(sum(pg[2:end] .* diff(x)), 1, rtol = 1e-2)
        @test isapprox(sum(pv[2:end] .* diff(x)), 1, rtol = 1e-2)
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

        @test isapprox(fwhm_l, fwhm(fl), rtol = 1e-2)
        @test isapprox(fwhm_g, fwhm(fg), rtol = 1e-2)
        @test isapprox(fwhm_v, LineBroadenings.fwhm_approx(fv), rtol = 1e-2) &&
              fwhm_v > fwhm(fl) &&
              fwhm_v > fwhm(fg)
    end
end
