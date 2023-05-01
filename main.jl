using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots
using Logging

Logging.disable_logging(Logging.Info)

include("hamiltonian.jl")
include("bond_convergence.jl")


plot_bond_convergence(4; spin=1, J=1.0, θ=0.0, maxiter=1000)

# Ψ = optimize_groundstate(N=Inf, spin=1, J=0.1, θ=0.0 , bond=4, maxiter=100)
# println(entanglement_spectrum(Ψ))
# entanglementplot(Ψ)

