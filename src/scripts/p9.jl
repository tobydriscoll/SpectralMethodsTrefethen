using CairoMakie, Printf, Polynomials, SpectralMethodsTrefethen
"p9 - polynomial interpolation in equispaced and Chebyshev pts"
function p9(N = 16)
    xx = range(-1, 1, 801)
    fig = Figure(size=(640, 300))
    cases = [
        ("equispaced points", [-1 + 2n / N for n in 0:N]),
        ("Chebyshev points", [cospi(n / N) for n in 0:N]),
    ]
    f(x) = 1 / (1 + 16x^2)
    for (i, (title, x)) in enumerate(cases)
        v = f.(x)
        ax = Axis(fig[1, i]; title, limits=(nothing, (-1, 1.5)))
        scatter!(ax, x, v)
        u = fit(Polynomial, x, v)       # interpolation
        uu = u.(xx)                     # evaluation of interpolant
        lines!(ax, xx, uu)
        maxerr = @sprintf("%.5g", maximum(abs, uu - f.(xx)))
        text!(ax, 0, -0.5; text="max error = $maxerr", align=(:center, :bottom))
    end
    return fig
end