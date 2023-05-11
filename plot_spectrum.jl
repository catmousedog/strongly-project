
"
Allow for plotting multiple entanglement spectra ontop of each other
"
function plot_entanglement(; spin=1, θ=0, bond=10, maxiter=500, clear=true,
    color=:blue, label=nothing, title="Entanglement Spectrum")

    ham = bilinear_biquadratic_hamiltonian(spin=spin, θ=θ)
    Ψ = optimize_groundstate(ham, bond=bond, maxiter=maxiter)
    spectrum = entanglement_spectrum(Ψ)
    if clear
        hline([1], color=:black, linewidth=1, linestyle=:dash, label=nothing, title=title)
    end
    scatter!(spectrum, yscale=:log10, xlabel="χ", ylabel="Schmidt coefficients",
        legend=!isnothing(label), label=label, xticks=1:length(spectrum), color=color)

    return spectrum
end

function annot(spectrum, χ, txt)
    annotate!(χ + 1, spectrum[χ], text(txt, :red, 10))
end

function linrange(start, stop, steps; endpoint=true)
    l = collect(range(start, stop=stop, length=steps))
    if !endpoint
        pop!(l)
    end
    return l
end

"
Plots entire spectrum for all θ in θ_range
"
function spectrum_approach(θ_range; spin=1, bond=12, maxiter=500, digits=4)
    steps = length(θ_range)
    spectra = Array{Float64}(undef, steps, bond)

    spectra[1, :] = plot_entanglement(;spin=spin, θ=θ_range[1], bond=bond, maxiter=maxiter, color=RGB(0., 0., 0.), 
    label="θ=$(round(θ_range[1], digits=digits))", title="θ=[$(round(θ_range[1], digits=digits)), $(round(last(θ_range), digits=digits)))")

    for (i, θ) in enumerate(θ_range[2:steps])
        c = i / (steps-1) * 0.9

        color = RGB(c, c, c)
        if i == steps-1
            label = "θ=$(round(last(θ_range), digits=digits))"
        else 
            label=nothing
        end

        spectra[i+1, :] = plot_entanglement(spin=spin, θ=θ, bond=bond, maxiter=maxiter, color=color, label=label, clear=false)
    end
    return spectra
end

"
Plot the change of the schmidt value χ as theta changes
"
function schmidt_range(θ_range, spectra, χ)
    scatter(θ_range, spectra[:, χ], yscale=:log10, xlabel="θ", ylabel="Schmidt coefficient $(χ)", legend=false, color=:black)
end