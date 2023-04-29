using Plots

include("hamiltonian.jl")

Ψ = optimize_groundstate(N=Inf, spin=1//2, J=1.0, θ=0 , max_bond=50, maxiter=100)
entanglementplot(Ψ)