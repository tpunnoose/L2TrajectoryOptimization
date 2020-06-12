# L2TrajectoryOptimization
EE364b Final Project

Tarun Punnoose and Nathan Kau

This is a ADMM-based trajectory optimization solver for problems with time-varying linear dynamics, quadratic state costs, control constraints, and group lasso (L2) cost on the control. The L2 control penalty encourages bang-off-bang control.

Please see "Paper.pdf" located in this repo for more detail.

## Using the code
To run the ADMM-based solver on an example satellite rendezvous problem:
```julia
include("src/Main.jl")
```

To run the CVX solver as a comparison (which doesn't converge):
```julia
include("src/CVXTest.jl")
```
