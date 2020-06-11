import Convex, ECOS

include("Dynamics.jl")
include("Parameters.jl")

##
function SelectState(t)
    return (6*(t-1)+1):(6*(t-1)+6)
end

function SelectControl(t)
    return (3*(t-1)+1):(3*(t-1)+3)
end

function DynamicsConstraints(X, U, parameters)::Vector{Convex.Constraint}
    constraints = []
    for t = 1:parameters.N
        push!(
            constraints,
            X[SelectState(t + 1)] ==
            CWDynamics(X[SelectState(t)], U[SelectState(t)], parameters),
        )
    end
    constraints
end

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

function CVX(parameters::Problem)
    N = parameters.N
    n = parameters.n
    m = parameters.m
    X = Convex.Variable(n * (N + 1))
    U = Convex.Variable(m * N)

    total_cost = L2Cost(U, parameters) + TerminalCost(X, parameters) + StageCost(X, parameters)
    constraints = DynamicsConstraints(X, U, parameters)
    prob = Convex.minimize(total_cost, constraints)
    X, U, prob
end
##
X, U, p = CVX(Problem(10))
Convex.solve!(p, ECOS.Optimizer(verbose=0))
