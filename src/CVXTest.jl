include("CVX.jl")
include("Parameters.jl")
include("ADMM.jl")

##
α = 50 # L2 weighting
ρ = 1e4
N = 100 # knot points

params = Problem(α, ρ, N)
X, U, p = CVX(params)
@time Convex.solve!(p, ECOS.Optimizer(verbose=1, maxit=200, feastol=1e-8))
Xs = StateTrajectoryToArray(X.value, params)
Us = ControlTrajectoryToArray(U.value, params)
U_norm = [norm(Us[:,i]) for i = 1:params.N]

println("\n\n -------------- CVX ------------")
@show p.optval
@show objective(X.value, U.value, params)

time = collect(1:N) .* params.Δt
##
a = plot(time, U_norm, label="Thrust")
# plot!(time, params.u_max * ones(N), label="control limit")
title!("L2 Group Lasso CVX")
xlabel!("Time (s)")
ylabel!("Thrust L2 Norm (N)")
png(a,"norm_cvx.png")
a
##
b = plot(time, Us', label=["u1" "u2" "u3"])
title!("L2 Group Lasso CVX")
xlabel!("Time (s)")
ylabel!("Thrust (N)")
png(b, "control_cvx.png")
b
##
c = plot(time, copy(Xs')[1:params.N,:], label=["px" "py" "pz" "dpx" "dpy" "dpz"])
title!("State Trajectory CVX")
xlabel!("Time (s)")
ylabel!("Position (m)")
png(c, "state_cvx.png")
c
