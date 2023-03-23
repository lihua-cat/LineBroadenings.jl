@doc raw"""

    Voigt(μ, σ, γ)

The *Voigt profile* is a probability distribution given by a convolution of a
Cauchy-Lorentz distribution and a Gaussian distribution.

# Definition
```math
f(x ;\mu, \sigma, \gamma) \equiv G(x;\mu, \sigma) \otimes L(x;\mu, \gamma)
```

It can be evaluated as

```math
f(x ;\mu, \sigma, \gamma)=\frac{\operatorname{Re}[w(z)]}{\sigma \sqrt{2 \pi}}
```

where $\mathrm{Re}[w(z)]$ is the real part of the *Faddeeva function* evaluated for

```math
z=\frac{x-\mu+i \gamma}{\sigma \sqrt{2}}
```

# Fields
- `μ::Real`: center
- `σ::Real`: σ in Gaussian distribution
- `γ::Real`: γ in Lorentz distribution

External links

* [Wikipedia](https://en.wikipedia.org/wiki/Voigt_profile)

"""
struct Voigt{T <: Real} <: Profile
    μ::T
    σ::T
    γ::T
    function Voigt{T}(µ::T, σ::T, γ::T) where {T <: Real}
        if σ >= zero(σ) && γ >= zero(γ)
            new{T}(µ, σ, γ)
        elseif σ < zero(σ)
            throw(DomainError(σ, "σ=$σ must be nonnegative"))
        else
            throw(DomainError(γ, "γ=$γ must be nonnegative"))
        end
    end
end

#### Outer constructors
Voigt(µ::T, σ::T, γ::T) where {T <: Real} = Voigt{T}(µ, σ, γ)
Voigt(μ::Real, σ::Real, γ::Real) = Voigt(promote(μ, σ, γ)...)
Voigt(μ::Integer, σ::Integer, γ::Integer) = Voigt(float(μ), float(σ), float(γ))

#### Parameters
params(d::Voigt) = (d.μ, d.σ, d.γ)

#### Statistics
w(z) = erfcx(-im * z) # Faddeeva function

function pdf(d::Voigt, x::Real)
    μ, σ, γ = params(d)
    if iszero(σ) && iszero(γ)
        p = x == μ ? Inf : zero(x)
    elseif iszero(σ)
        dl = Lorentz(μ, γ)
        p = pdf(dl, x)
    elseif iszero(γ)
        dg = Gaussian(μ, σ)
        p = pdf(dg, x)
    else
        z = (x - μ + im * γ) / (√(2) * σ)
        p = real(w(z)) / (SQRT2π * σ)
    end
    return p
end

peak(d::Voigt) = pdf(d, d.μ)

function fwhm_approx(d::Voigt{T}) where {T}
    fl = W_L * d.γ
    fg = W_G * d.σ
    return 0.5346fl + √(0.2166fl^2 + fg^2) |> T
end
