using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots

include("hamiltonian.jl")
include("groundstate.jl")
include("plot_spectrum.jl")

# using Logging
# disable_logging(LogLevel(Info))


H = HAFM_staggered(spin=1, J=1.0, g=0.01)
Ψ_HAFM = optimize_groundstate(H, bond=12, maxiter=500)
spectrum_HAFM = entanglement_spectrum(Ψ_HAFM)
plot_entanglement(spectrum_HAFM)
display(current())
