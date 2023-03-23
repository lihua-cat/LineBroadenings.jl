@doc raw"""

    Lorentz(μ, γ)

Lorentz distribution, also called Cauchy distribution.

```math
f\left(x ; \mu, \gamma\right)=\frac{1}{\pi \gamma\left[1+\left(\frac{x-\mu}
{\gamma}\right)^{2}\right]}=\frac{1}{\pi \gamma}\left[\frac{\gamma^{2}}
{\left(x-\mu\right)^{2}+\gamma^{2}}\right]
```

# Fields
- `μ::Real`: location parameter
- `γ::Real`: scale parameter

External links

* [Wikipedia](https://en.wikipedia.org/wiki/Cauchy_distribution)

"""
struct Lorentz{T <: Real} <: Profile
    μ::T
    γ::T
    function Lorentz{T}(µ::T, γ::T) where {T <: Real}
        γ >= zero(γ) ? new{T}(µ, γ) : throw(DomainError(γ, "γ=$γ must be nonnegative"))
    end
end

#### Outer constructors
Lorentz(µ::T, γ::T) where {T <: Real} = Lorentz{T}(µ, γ)
Lorentz(μ::Real, γ::Real) = Lorentz(promote(μ, γ)...)
Lorentz(μ::Integer, σ::Integer) = Lorentz(float(μ), float(σ))

#### Parameters
params(d::Lorentz) = (d.μ, d.γ)

#### Statistics
const W_L = 2

function pdf(d::Lorentz, x::Real)
    μ, γ = params(d)
    if iszero(γ)
        p = x == μ ? Inf : zero(x)
    else
        p = 1 / π / γ * (γ^2 / ((x - μ)^2 + γ^2))
    end
    return p
end

peak(d::Lorentz) = 1 / π / d.γ

fwhm(d::Lorentz) = W_L * d.γ
