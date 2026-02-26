using CairoMakie, LaTeXStrings
using SpectralMethodsTrefethen, SpecialFunctions
"""
p30 - spectral integration, ODE style (compare p12)\n
   use `chebweights`, `clencurt`, or `gauss` as the argument
"""
function p30(quadrule = chebweights)
    # Computation: various values of N, four functions:
    Nmax = 50
    funs = [x -> abs(x)^3    x -> exp(-x^(-2));
            x -> 1 / (1 + x^2)     x -> x^10]
    ints = [1/2 2*(exp(-1) + sqrt(π) * (erf(1) - 1));   π/2 2/11]
    E = zeros(2, 2, Nmax)
    for N = 1:Nmax
        x, w = quadrule(N)
        for i in CartesianIndices(funs)
            f = funs[i].(x)
            E[i, N] = abs(dot(w, f) - ints[i])
        end
    end

    # Plot results:
    fig = Figure()
    titles = [L"|x|^3"  L"\exp(-x^2)"; L"1/(1+x^2)"  L"x^{10}"]
    ax = [Axis(fig[i, j]; title=titles[i, j], yscale=log10) for i in 1:2, j in 1:2]
    linkaxes!(ax...)
    @. E[iszero(E)] = NaN
    for i in 1:2, j in 1:2
        scatterlines!(ax[i, j], 1:Nmax, E[i, j, :]; markersize=6)
    end
    ax[1, 1].ylabel = ax[2, 1].ylabel = "error"
    ax[2, 1].xlabel = ax[2, 2].xlabel = L"N"
    return fig
end

function chebweights(N)
    i = 1:N
    D, x = cheb(N)
    x = x[i]
    w = inv(D[i, i])[1, :]
    return x, w
end

p30()
