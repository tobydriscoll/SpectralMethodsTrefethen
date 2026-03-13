using CairoMakie, Polynomials
"p10 - polynomials and corresponding equipotential curves"
function p10(N=16)
    fig = Figure()
    cases = (
        ("equispaced points", [-1 + 2n / N for n in 0:N]),
        ("Chebyshev points", [cospi(n / N) for n in 0:N]),
    )
    for (i, (title, x)) in enumerate(cases)
        p = fromroots(x)

        # Plot p(x) over [-1,1]:
        scatter(fig[i, 1], x, zeros(size(x)), axis=(; title))
        lines!(fig[i, 1], -1..1, x -> p(x))

        # Plot equipotential curves:
        xx, yy = -1.3:0.02:1.3, -1:0.02:1
        P = [p(x + im * y) for x in xx, y in yy]
        logP = @. log10(max(abs(P), 1e-16))
        ax = Axis(fig[i, 2]; title, aspect=DataAspect())
        contour!(ax, xx, yy, logP; levels=-4:0, linewidth=2,
            colormap=:viridis, colorrange=(-4, 0))
        scatter!(ax, x, zeros(size(x)); color=:black)
    end
    return fig
end