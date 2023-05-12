using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots

include("hamiltonian.jl")
include("groundstate.jl")
include("plot_spectrum.jl")

using Logging
disable_logging(LogLevel(Info))


# ham_p = bilinear_biquadratic_hamiltonian_perturbed(spin=1, θ=0, g=0.2)
# Ψ_p = optimize_groundstate(ham_p, bond=12, maxiter=500)
# entanglementplot(Ψ_p)
# display(current())
