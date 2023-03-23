@doc raw"""

    Gaussian(μ, σ)

Gaussian distribution, also called normal distribution.

# Definition
```math
f\left(x ; \mu, \sigma\right)=\frac{1}{\sigma \sqrt{2 \pi}} \exp{\left[-\frac{1}{2}
\left(\frac{x-\mu}{\sigma}\right)^{2}\right]}
```

# Fields
- `μ::Real`: mean or expectation
- `σ::Real`: standard deviation

External links

* [Wikipedia](https://en.wikipedia.org/wiki/Normal_distribution)

"""
struct Gaussian{T <: Real} <: Profile
    μ::T
    σ::T
    function Gaussian{T}(µ::T, σ::T) where {T <: Real}
        σ >= zero(σ) ? new{T}(µ, σ) : throw(DomainError(σ, "σ=$σ must be nonnegative"))
    end
end

#### Outer constructors
Gaussian(µ::T, σ::T) where {T <: Real} = Gaussian{T}(µ, σ)
Gaussian(μ::Real, σ::Real) = Gaussian(promote(μ, σ)...)
Gaussian(μ::Integer, σ::Integer) = Gaussian(float(μ), float(σ))

#### Parameters
params(d::Gaussian) = (d.μ, d.σ)

#### Statistics
const W_G = 2 * √(2 * log(2))
const SQRT2π = √(2π)

function pdf(d::Gaussian, x::Real)
    μ, σ = params(d)
    if iszero(σ)
        p = x == μ ? Inf : zero(x)
    else
        p = 1 / σ / SQRT2π * exp(-((x - μ) / σ)^2 / 2)
    end
    return p
end

peak(d::Gaussian) = 1 / d.σ / SQRT2π

fwhm(d::Gaussian) = W_G * d.σ
