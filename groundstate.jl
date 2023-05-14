
function optimize_groundstate(ham; bond::Int64=12, maxiter::Int64=100)
    Ψ = InfiniteMPS(ℂ^Int(2 * ham.spin + 1), ℂ^(bond))
    Ψ, = find_groundstate(Ψ, ham.H, VUMPS(maxiter=maxiter))
    return Ψ
end

function linrange(start, stop, steps; endpoint=true)
    l = collect(range(start, stop=stop, length=steps))
    if !endpoint
        pop!(l)
    end
    return l
end

"
Give the ground states of the bilinear_biquadratic_hamiltonian for each θ in θ_range
"
function bilinear_biquadratic_θ_range(θ_range; spin=1, J=1.0, bond::Int64=12, maxiter::Int64=100)
    return [
        begin
            ham = bilinear_biquadratic_hamiltonian(spin=spin, J=J, θ=θ)
            optimize_groundstate(ham, bond=bond, maxiter=maxiter)
        end for θ in θ_range
    ]
end

"
Calculate spectrum for multiple InfiniteMPS
"
function entanglement_spectra(Ψ_list)
    spectra = nothing
    for (i, Ψ) in enumerate(Ψ_list)
        if i == 1 # don't know how to get bond dimension from InfiniteMPS so have to do this
            spectrum = entanglement_spectrum(Ψ)
            spectra = zeros(length(Ψ_list), length(spectrum))
            spectra[1, :] = spectrum
        else
            spectra[i, :] = entanglement_spectrum(Ψ)
        end
    end
    return spectra
end

function correlation_lengths(Ψ_list)
    return [correlation_length(Ψ) for Ψ in Ψ_list]
end

