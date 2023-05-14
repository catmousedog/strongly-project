
"
Allow for plotting multiple entanglement spectra ontop of each other
"
function plot_entanglement(spectrum; clear=true, color=:blue, label=nothing, title="Entanglement Spectrum")
    # to clear plot we draw hline at 1
    if clear
        hline([1], color=:black, linewidth=1, linestyle=:dash, label=nothing, title=title)
    end

    scatter!(spectrum, yscale=:log10, xlabel="χ", ylabel="Schmidt coefficients",
        legend=true, label=label, xticks=1:length(spectrum), color=color)

    return spectrum
end

"
Plots all spectra uses θ_range for labels and title
"
function plot_spectra(θ_range, spectra; digits=4)
    @assert length(θ_range) == length(Ψ_list)

    steps = length(θ_range)

    for i in 1:steps
        θ = θ_range[i]
        spectrum = spectra[i, :]

        c = i / (steps - 1) * 0.9
        color = RGB(c, c, c)

        if i == 1
            plot_entanglement(spectrum, color=color, label="θ=$(round(θ, digits=digits))", clear=true)
        elseif i == steps
            plot_entanglement(spectrum, color=color, label="θ=$(round(θ, digits=digits))", clear=false)
        else
            plot_entanglement(spectrum, color=color, clear=false)
        end
    end

    title!("θ=[$(round(θ_range[1], digits=digits)), $(round(last(θ_range), digits=digits))]")
end

"
Plot the change of the schmidt value χ as theta changes
"
function plot_schmidt_range(θ_range, spectra, χ)
    scatter(θ_range, spectra[:, χ], yscale=:log10, xlabel="θ", ylabel="Schmidt coefficient $(χ)", legend=false, color=:black)
end

function plot_correlations(θ_range, correlations)
    scatter(θ_range, correlations, yscale=:log10, xlabel="θ", ylabel="correlation lengths", legend=false, color=:black)
end