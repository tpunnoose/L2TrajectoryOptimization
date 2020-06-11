using LinearAlgebra
using Plots

include("Parameters.jl")
include("ADMM.jl")
include("Utils.jl")

ρ = 1.0
α = 1e-3

p = Problem(α, ρ)

(X,U) = ADMM(p, 10)

group_norm = [norm(U[Select3(i)]) for i=1:p.N]

plot(1:p.N, group_norm)
plot!(xlabel="Time", ylabel="Thrust Norm", title="L2 Group Norm ADMM")
