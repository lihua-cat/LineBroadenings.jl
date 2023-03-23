abstract type Profile end

Broadcast.broadcastable(d::Profile) = Ref(d)

include("gaussian.jl")
include("lorentz.jl")
include("voigt.jl")
