include("CVX.jl")
include("Parameters.jl")
include("ADMM.jl")

##
α = 1 # L2 weighting
ρ = 1
N = 100 # knot points

params = Problem(α, ρ, N)
X, U, p = CVX(params)
@time Convex.solve!(p, ECOS.Optimizer(verbose=0))
Xs = StateTrajectoryToArray(X.value, params)
Us = ControlTrajectoryToArray(U.value, params)
U_norm = [norm(Us[:,i]) for i = 1:params.N]

println("\n\n -------------- CVX ------------")
@show p.optval
@show objective(X.value, U.value, params)


##
plot(U_norm)
title!("CVX Control norm")
##
plot(Us')
title!("CVX Control")

##
plot(copy(Xs'))
title!("CVX State Trajectory")
