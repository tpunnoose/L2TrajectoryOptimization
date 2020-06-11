using LinearAlgebra
using Plots

include("Parameters.jl")
include("ADMM.jl")
include("Utils.jl")

ρ = 1
α = 1

p = Problem(α, ρ, 100)

X, U, Y = ADMM(p, 100)

u_group_norm = [norm(U[SelectControl(i)]) for i = 1:p.N]
y_group_norm = [norm(Y[SelectControl(i)]) for i = 1:p.N]

##
@show objective(X, U, p)

##
a = plot(1:p.N, u_group_norm)
plot!(1:p.N, y_group_norm)
plot!(xlabel = "Time", ylabel = "Thrust Norm", title = "L2 Group Norm ADMM")

##
b = plot(reshape(X, p.n, p.N + 1)')
title!("ADMM State Trajectory")

# objective()
