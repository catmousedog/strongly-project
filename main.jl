using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots


include("hamiltonian.jl")

Ψ = optimize_groundstate(N=Inf, spin=1, J=0.1, θ=0.0 , bond=4, maxiter=100)
entanglementplot(Ψ)

