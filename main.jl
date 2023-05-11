using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots

include("hamiltonian.jl")
include("plot_spectrum.jl")

using Logging
disable_logging(LogLevel(Info))



θ_range = linrange(0, atan(1/3), 10; endpoint=false)

spectra = spectrum_approach(θ_range, spin=1, bond=12)

schmidt_range(θ_range, spectra, 4)
