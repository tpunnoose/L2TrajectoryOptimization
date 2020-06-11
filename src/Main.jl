using LinearAlgebra
using Plots

include("Parameters.jl")
include("ADMM.jl")
include("CVX.jl")
include("Utils.jl")

ρ = 100
α = 5

p = Problem(α, ρ, 100)

X, U, Y, Z = ADMMGroupConsensus(p, 100)

u_group_norm = [norm(U[SelectControl(i)]) for i = 1:p.N]
y_group_norm = [norm(Y[SelectControl(i)]) for i = 1:p.N]
z_group_norm = [norm(Z[SelectControl(i)]) for i = 1:p.N]

##
@show objective(X, U, p)

##
a = plot(1:p.N, u_group_norm, label="LQR")
plot!(1:p.N, y_group_norm, label="L2")
plot!(1:p.N, z_group_norm, label="Control Constraint")
plot!(xlabel = "Time", ylabel = "Thrust Norm", title = "L2 Group Norm ADMM")

##
b = plot(reshape(X, p.n, p.N + 1)')
title!("ADMM State Trajectory")

# objective()
