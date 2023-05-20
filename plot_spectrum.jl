
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
Plots all spectra uses range for labels and title
"
function plot_spectra(range, range_name, spectra; digits=4)
    steps = length(range)

    for i in 1:steps
        x = range[i]
        spectrum = spectra[i, :]

        c = i / steps * 0.9
        color = RGB(c, c, c)

        if i == 1
            plot_entanglement(spectrum, color=color, label="$(range_name)=$(round(x, digits=digits))", clear=true)
        elseif i == steps
            plot_entanglement(spectrum, color=color, label="$(range_name)=$(round(x, digits=digits))", clear=false)
        else
            plot_entanglement(spectrum, color=color, clear=false)
        end
    end

    title!("$(range_name)=[$(round(range[1], digits=digits)), $(round(last(range), digits=digits))]")
end
