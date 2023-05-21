
function optimize_groundstate(ham; bond::Int64=12, maxiter::Int64=100)
    data = fill(TensorMap(rand,ComplexF64,ℂ^bond*ℂ^Int(2 * ham.spin + 1),ℂ^bond), ham.unit_cells);
    Ψ = InfiniteMPS(data);
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

function logticks(sequence)
    min_val = minimum(sequence)
    max_val = maximum(sequence)
    min_log = floor(log10(min_val))
    max_log = ceil(log10(max_val))
    ticks = [10.0^i for i in min_log:max_log]
    return ticks
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
Give the ground state for a given hamiltonian 'ham' for each bond in bond_range
"
function correlation_bond_range(bond_range, ham; maxiter::Int64=500)
    return [
        begin
            Ψ = optimize_groundstate(ham, bond=bond, maxiter=maxiter)
            correlation_length(Ψ)
        end for bond in bond_range
    ]
end


"
Give the ground state of the petrurbed hamiltonian H for each g in g_range
"
function perturbation_range(H_perturbed, g_range; spin=1, J=1.0, bond::Int64=12, maxiter::Int64=100)
    return [
        begin
            ham = H_perturbed(spin=spin, J=J, g=g)
            optimize_groundstate(ham, bond=bond, maxiter=maxiter)
        end for g in g_range
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
