using LinearAlgebra

include("Parameters.jl")
include("ADMM.jl")

p = Problem()

ρ = 1.0

Y = ones(p.N*p.m)
Λ = 0.1*ones(p.N*p.m)

ADMM(p, 10)
