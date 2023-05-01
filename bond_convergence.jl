
function last_singular(;N=Inf, spin=1, J=1.0, θ=0.0 , bond=30, maxiter=600)
    Ψ = optimize_groundstate(N=N, spin=spin, J=J, θ=θ , bond=bond, maxiter=maxiter)
    # lowest singular value determines how well the chosen bond dimension is
    eigv = entanglement_spectrum(Ψ)
    return eigv[end]
end

function plot_bond_convergence(max_bond; spin=1, J=1.0, θ=0.0, maxiter=600)
    bonds = 2:max_bond

    # plot last singular values for each bond dimension
    last_spectrum = [last_singular(spin=spin, J=J, θ=θ, bond=bond, maxiter=maxiter) for bond in bonds]
    scatter(bonds, last_spectrum, yscale=:log10, label="last singular value")

    # plot entire spectrum up for max bond dimension
    Ψ = optimize_groundstate(N=Inf, spin=spin, J=J, θ=θ, bond=max_bond, maxiter=maxiter)
    full_spectrum = entanglement_spectrum(Ψ)
    full_spectrum = full_spectrum[2:end]
    scatter!(bonds, full_spectrum, yscale=:log10, label="full spectrum")

    xticks!(Int.(bonds))
    xlabel!("Bond Dimension")
    ylabel!("Entanglement Spectrum")
    title!("spin=$(spin),  J=$(J),  θ=$(θ),  max bond=$(max_bond)")
end
