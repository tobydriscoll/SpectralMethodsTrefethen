using CairoMakie, LaTeXStrings
using ForwardDiff, SpectralMethodsTrefethen
"p7 - accuracy of periodic spectral differentiation"
function p7(Nmax=50)
    # Compute derivatives for various values of N:
    N = 6:2:Nmax
    funs = [x -> abs(sin(x))^3            x -> exp(-sin(x / 2)^(-2));
            x -> 1 / (1 + sin(x / 2)^2)   x -> sin(10x)]
    E = zeros(2, 2, length(N))
    for (k, N) in enumerate(N)
        x, D, _ = fourier(N)
        for (i, f) in pairs(funs)
            v = f.(x)
            vprime = ForwardDiff.derivative.(f, x)
            E[i, k] = maximum(abs, D * v - vprime)
        end
    end

    # Plot results:
    titles = [L"|\sin(x)|^3" L"\exp(-\sin^{-2}(x/2))"; L"1/(1+\sin^2(x/2))" L"\sin(10x)"]
    fig = Figure()
    ax = [Axis(fig[i, j]; title=titles[i, j], yscale=log10) for i in 1:2, j in 1:2]
    linkaxes!(ax...)
    for (i, ax) in pairs(ax)
        scatterlines!(ax, N, E[i, :])
    end
    ax[1, 1].ylabel = ax[2, 1].ylabel = "error"
    ax[2, 1].xlabel = ax[2, 2].xlabel = L"N"
    return fig
end