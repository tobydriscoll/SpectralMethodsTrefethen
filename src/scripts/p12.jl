using CairoMakie, LaTeXStrings, SpectralMethodsTrefethen
"p12 - accuracy of Chebyshev spectral differentiation (compare p7)"
function p12(Nmax = 50)
    # Compute derivatives for various values of N:
    # 3rd deriv in BV, C^∞; analytic, polynomial
    funs = [x -> abs(x)^3            x -> exp(-x^(-2));
            x -> 1 / (1 + x^2)       x -> x^10]
    E = zeros(2, 2, Nmax)
    for N in 1:Nmax
        D, x = cheb(N)
        for (i, f) in pairs(funs)
            v = f.(x)
            df_dx = ForwardDiff.derivative.(f, x)
            E[i, N] = maximum(abs, D * v - df_dx)
        end
    end

    # Plot results:
    fig = Figure()
    titles = [L"|x|^3"  L"\exp(-x^2)"; L"1/(1+x^2)"  L"x^{10}"]
    ax = [Axis(fig[i, j]; title=titles[i, j], yscale=log10) for i in 1:2, j in 1:2]
    linkaxes!(ax...)
    for (i, ax) in pairs(ax)
        scatterlines!(ax, 1:Nmax, E[i, :])
    end
    ax[1, 1].ylabel = ax[2, 1].ylabel = "error"
    ax[2, 1].xlabel = ax[2, 2].xlabel = L"N"
    return fig
end