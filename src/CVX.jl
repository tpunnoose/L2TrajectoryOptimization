import Convex, ECOS

include("Parameters.jl")

##
function SelectState(t)
    return (6*(t-1)+1):(6*(t-1)+6)
end

function SelectControl(t)
    return (3*(t-1)+1):(3*(t-1)+3)
end

function DynamicsConstraints(X, U, parameters::Problem)::Vector{Convex.Constraint}
    constraints = []
    for t = 1:parameters.N
        push!(
            constraints,
            X[SelectState(t + 1)] == parameters.A * X[SelectState(t)] + parameters.B * U[SelectControl(t)]
        )
    end
    constraints
end

BoxConstraint(U, parameters) = [abs(U) <= 1e-2]

InitialConstraint(X, parameters::Problem) = [X[SelectState(1)] == parameters.x0]

function TerminalCost(X, parameters::Problem)
    final_x = X[SelectState(parameters.N + 1)]
    cost = Convex.quadform(final_x, parameters.Q_f)
end

function StageCost(X, parameters::Problem)
    cost = 0
    for t = 1:parameters.N
        cost += Convex.quadform(X[SelectState(t)], parameters.Q_k)
    end
    cost
end

function L2Cost(U, parameters::Problem)
    cost = 0
    for t = 1:parameters.N
        cost += Convex.norm(U[SelectControl(t)], 2)
    end
    cost
end

function TotalCost(X, U, parameters)
    parameters.α * L2Cost(U, parameters) + StageCost(X, parameters) + TerminalCost(X, parameters)
end

function CVX(parameters::Problem, alpha=0.1)
    N = parameters.N
    n = parameters.n
    m = parameters.m
    X = Convex.Variable(n * (N + 1))
    U = Convex.Variable(m * N)

    total_cost = TotalCost(X, U, parameters)
    constraints = [DynamicsConstraints(X, U, parameters); InitialConstraint(X, parameters); BoxConstraint(U, parameters)]
    prob = Convex.minimize(total_cost, constraints)
    X, U, prob
end
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
