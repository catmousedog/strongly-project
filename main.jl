using MPSKitModels, MPSKit, TensorKit, TensorOperations, Plots

include("hamiltonian.jl")
include("groundstate.jl")
include("plot_spectrum.jl")

# using Logging
# disable_logging(LogLevel(Info))


"
Don't keep MPS in memory
"
θ_range = linrange(0, π/4, 20, endpoint=false)
corr = []
for θ in θ_range
    H = bilinear_biquadratic_hamiltonian(spin=1, θ=θ)
    Ψ = optimize_groundstate(H, bond=24, maxiter=500) # good bond dimension for convergence
    push!(corr, correlation_length(Ψ))
end
plot_correlations(θ_range, corr)


"
fixed θ
varying bond
"
# bonds = 12:24
# correlations = []
# for D in bonds
#     println(D)
#     H = bilinear_biquadratic_hamiltonian(spin=1, θ=π/4)
#     Ψ = optimize_groundstate(H, bond=D, maxiter=500)
#     push!(correlations, correlation_length(Ψ))
# end
# scatter(bonds, correlations, yscale=:log10, xlabel="bond", ylabel="correlation lengths", legend=false, color=:black)


"
Keep MPS in memory
"
# θ_Nematic_approach = linrange(0, π/4, 20; endpoint=false)
# Ψ_Nematic_list = bilinear_biquadratic_θ_range(θ_Nematic_approach, spin=1, bond=20, maxiter=500)
# correlations_Nematic = correlation_lengths(Ψ_Nematic_list)
# plot_correlations(θ_Nematic_approach, correlations_Nematic)

display(current())