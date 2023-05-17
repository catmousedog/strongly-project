using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots

include("hamiltonian.jl")
include("groundstate.jl")
include("plot_spectrum.jl")

# using Logging
# disable_logging(LogLevel(Info))


HAFM = HAFM_staggered()
Ψ_HAFM = optimize_groundstate(HAFM, bond=12, maxiter=500)
spectrum_HAFM = entanglement_spectrum(Ψ_HAFM)
plot_entanglement(spectrum_HAFM)
display(current())
