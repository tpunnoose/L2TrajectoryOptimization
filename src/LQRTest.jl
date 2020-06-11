import Convex, ECOS
using Plots

include("CVX.jl")
include("ADMM.jl")
##
α = 1e-2
ρ = 1
N = 100
p = Problem(α, ρ, N)

A = ones(p.m * p.N) + randn(p.m * p.N)*1e-4
Y = ones(p.m * p.N) + randn(p.m * p.N)*1e-4

LQRCVX(p.ρ, Y, A, p)
Xs = StateTrajectoryToArray(X.value, p)
Us = ControlTrajectoryToArray(U.value, p)
a1 = plot(Xs')
title!(a1, "CVX State")
b1 = plot(Us')
title!(b1, "CVX Control")
plot(a1, b1)

X2, U2 = LQR(p.ρ, Y, A, p)
Xs2 = StateTrajectoryToArray(X2, p)
Us2 = ControlTrajectoryToArray(U2, p)
a2 = plot(Xs2')
title!(a2, "LQR State")
b2 = plot(Us2')
title!(b2, "LQR Control")
plot(a2, b2)

##
plot(Xs'[:,1])
plot!(Xs2'[:,1])
# plot(Xs')
# plot!(Xs2')
title!("LQR vs CVX State")

##
plot(Us'[:,3])
plot!(Us2'[:,3])
title!("LQR vs CVX Control")
