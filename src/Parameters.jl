using LinearAlgebra

mutable struct Problem
    N::Int64 # number of nodes
    n::Int64 # number of states
    m::Int64 # number of controls
    x0::Vector{Float64} # initial state
    xf::Vector{Float64} # final state
    Δt::Float64
    Q_f::Matrix{Float64} # terminal state penalty
    Q_k::Matrix{Float64} # stagewise cost penalty
    n_::Float64 # mean motion
    m_ego::Float64 # mass of controlled spacecraft
    A::Matrix{Float64} # dynamics state matrix
    B::Matrix{Float64} # control state matrix
    α::Float64 # regularizer parameter
    ρ::Float64 # augmented lagrangian parameter
    u_max::Float64 # maximum control magnitude
end

function Problem(α, ρ, N=100)
    μ = 3.99 * 10^14 # Standard gravitational parameter m^3 s^-2
    a = 6731.0 * 10^3 # Earth radius in m
    alt = 5e5 # Altitude of the target satellite in m
    orbit_radius = a + alt

    n_ = sqrt(μ / orbit_radius^3)

    # N = 100 # Number of node points
    n = 6 # State vector dimension
    m = 3 # Control vector dimension
    tf = 6000 # Trajectory duration in s (approx. = 1 revolutions)
    Δt = tf / (N - 1)

    Q_f = 1e5 * Matrix{Float64}(I, n, n)
    Q_k = 1e-3* Matrix{Float64}(I, n, n)

    Δr0 = [20.0, 30.0, 10.0] # initial position delta
    Δrd0 = [0.4, 0.6, -0.1] # intial velocity delta
    # Δr0 = [0.0, 0.0, 0.0] # initial position delta
    # Δrd0 = [0.0, 0.0, 0.0] # intial velocity delta
    x0 = [Δr0..., Δrd0...] # initial state
    xf = zeros(n)
    m_ego = 4.6

    A_c = zeros(n,n)
    A_c[1:3, 4:6] = Matrix{Float64}(I, 3, 3)
    A_c[4, 5] = 2*n_
    A_c[5, 2] = 3*n_^2
    A_c[5, 4] = -2*n_
    A_c[6, 3] = -n_^2

    B_c = zeros(n,m)
    B_c[4:6, :] = Matrix{Float64}(I, 3, 3)/m_ego

    cont_sys = zeros(n+m, n+m)
    cont_sys[1:n, 1:n] = A_c
    cont_sys[1:n, (n+1):(n+m)] = B_c

    disc_sys = exp(cont_sys*Δt)
    A_d = disc_sys[1:n, 1:n]
    B_d = disc_sys[1:n, (n+1):(n+m)]

    u_max = 100

    return Problem(N, n, m, x0, xf, Δt, Q_f, Q_k, n_, m_ego, A_d, B_d, α, ρ, u_max)
end
