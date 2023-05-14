using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots

include("hamiltonian.jl")
include("groundstate.jl")
include("plot_spectrum.jl")

using Logging
disable_logging(LogLevel(Info))

H = bilinear_biquadratic_hamiltonian_perturbed(spin=1, θ=0.0, J=1.0, g=0.1)
Ψ = optimize_groundstate(H, bond=12, maxiter=100)
spectrum = entanglement_spectrum(Ψ)
plot_entanglement(spectrum)

# θ_Nematic_approach = linrange(0, π/4, 20; endpoint=false)
# Ψ_Nematic_list = bilinear_biquadratic_θ_range(θ_Nematic_approach, spin=1, bond=12, maxiter=500)
# correlations_Nematic = correlation_lengths(Ψ_Nematic_list)
# plot_correlations(θ_Nematic_approach, correlations_Nematic)

display(current())
