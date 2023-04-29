using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots

include("hamiltonian.jl")

Ψ = optimize_groundstate(N=Inf, spin=1, J=0.1, θ=π , max_bond=30, maxiter=100)
entanglementplot(Ψ)