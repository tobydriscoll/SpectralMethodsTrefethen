using CairoMakie, SpectralMethodsTrefethen, SpecialFunctions
"p31 - gamma function via complex integral, trapezoid rule"
function p31(N = 70, c = -11, r = 16)
    # points on circle of integration:
    t = [c + r * cispi(q) for q in 2 * (0.5:N-0.5) / N]
    # trapezoid integration:
    gaminv(z) = sum(exp(t) * t^(-z) * (t - c) for t in t) / N

    x, y = range(-3.5, 4, 121), range(-2.5, 2.5, 121)
    Γ = [1 / gaminv(complex(x, y)) for x in x, y in y]
    fig = Figure(size=(700, 500))
    ax = Axis3(fig[1, 1]; xlabel="Re(z)", ylabel="Im(z)", zlabel="|Γ(z)|")
    G = @. min(6, abs(Γ))
    plt = surface!(ax, x, y, G; colormap=:lajolla, colorscale=log2, colorrange=(1e-2, 6))
    Colorbar(fig[1, 2], plt)
    return fig
end

p31()
