using LinearAlgebra
using Plots

include("Parameters.jl")
include("ADMM.jl")
include("CVX.jl")
include("Utils.jl")

ρ = 20000
α = 50

p = Problem(α, ρ, 100)

X, U, Y, Z = ADMMGroupConsensus(p, 1000)

u_group_norm = [norm(U[SelectControl(i)]) for i = 1:p.N]
y_group_norm = [norm(Y[SelectControl(i)]) for i = 1:p.N]
z_group_norm = [norm(Z[SelectControl(i)]) for i = 1:p.N]
time = collect(1:p.N) .* params.Δt

##
@show objective(X, U, p)

##
a = plot(time, u_group_norm, label="Thrust")
# scatter!(a, time, u_group_norm, label="LQR")
# plot!(time, y_group_norm, label="L2")
# plot!(time, z_group_norm, label="Control Constraint")
plot!(xlabel = "Time (s)", ylabel = "Thrust L2 Norm (N)", title = "L2 Group Lasso ADMM")
png(a,"norm_admm.png")
a
##
Us = ControlTrajectoryToArray(U, p)
b = plot(time, Us',label=["u1" "u2" "u3"])
title!("L2 Group Lasso ADMM")
xlabel!("Time")
ylabel!("Thrust (N)")
png(b,"control_admm.png")
b

##
c = plot(time, reshape(X, p.n, p.N + 1)'[1:p.N,:], label=["px" "py" "pz" "dpx" "dpy" "dpz"])
title!("State Trajectory ADMM")
xlabel!("Time (s)")
ylabel!("Position (m)")
png(c,"state_admm.png")
c
