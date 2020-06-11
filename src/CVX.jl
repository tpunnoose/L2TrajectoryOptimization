import Convex, ECOS

include("Parameters.jl")
include("Utils.jl")

##
function SelectState(t)
    return (6*(t-1)+1):(6*(t-1)+6)
end

function SelectControl(t)
    return (3*(t-1)+1):(3*(t-1)+3)
end

function DynamicsConstraints(X, U, parameters::Problem)
    constraints = Vector{Convex.Constraint}()
    for t = 1:parameters.N
        push!(
            constraints,
            X[SelectState(t + 1)] ==
            parameters.A * X[SelectState(t)] +
            parameters.B * U[SelectControl(t)],
        )
    end
    constraints
end

BoxConstraint(U, parameters) = [abs(U) <= parameters.u_max]

InitialConstraint(X, parameters::Problem) = [X[SelectState(1)] == parameters.x0]

function TerminalCost(X, parameters::Problem)
    final_x = X[SelectState(parameters.N + 1)]
    cost = 0.5 * Convex.quadform(final_x, parameters.Q_f)
end

function StageCost(X, parameters::Problem)
    cost = 0
    for t = 1:parameters.N
        cost += 0.5 * Convex.quadform(X[SelectState(t)], parameters.Q_k)
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
    parameters.α * L2Cost(U, parameters) +
    StageCost(X, parameters) +
    TerminalCost(X, parameters)
end

function CVX(parameters::Problem, alpha = 0.1)
    N = parameters.N
    n = parameters.n
    m = parameters.m
    X = Convex.Variable(n * (N + 1))
    U = Convex.Variable(m * N)

    total_cost = TotalCost(X, U, parameters)
    constraints = [
        DynamicsConstraints(X, U, parameters)
        InitialConstraint(X, parameters)
        BoxConstraint(U, parameters)
    ]
    prob = Convex.minimize(total_cost, constraints)
    X, U, prob
end


function ControlCost(U, Y, Λ̄, params)
    R = params.ρ / 2 * I(params.m)
    U_ref = Y - Λ̄
    cost = 0
    for t = 1:params.N
        e = U[SelectControl(t)] - U_ref[SelectControl(t)]
        cost += 0.5 * Convex.quadform(e, R)
    end
    return cost
end

function LQRCVX(ρ, Y, Λ̄, p)
    X, U = Convex.Variable(p.n * (p.N + 1)), Convex.Variable(p.m * p.N)
    cost = ControlCost(U, Y, Λ̄, p) + StageCost(X, p) + TerminalCost(X, p)
    constraints = [DynamicsConstraints(X, U, p); InitialConstraint(X, p)]

    problem = Convex.minimize(cost, constraints)
    Convex.solve!(problem, ECOS.Optimizer())
    return X.value, U.value
end
