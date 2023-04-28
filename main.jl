include("hamiltonian.jl")

using .hamiltonian
using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots

Ψ = optimize_groundstate(N=Inf, spin=1, J=1.0, θ=0.0, max_bond=50, maxiter=100)
entanglementplot(Ψ)