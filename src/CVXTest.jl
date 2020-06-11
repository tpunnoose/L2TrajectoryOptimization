include("CVX.jl")
include("Parameters.jl")

##
α = 0.001 # L2 weighting
ρ = 1
N = 100 # knot points

params = Problem(α, ρ, N)
X, U, p = CVX(params)
Convex.solve!(p, ECOS.Optimizer(verbose=0))
Xs = StateTrajectoryToArray(X.value, params)
Us = ControlTrajectoryToArray(U.value, params)
U_norm = [norm(Us[:,i]) for i = 1:params.N]

##
plot(U_norm)
title!("Control norm")

##
plot(copy(Xs'))
title!("State Trajectory")
