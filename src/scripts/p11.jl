using CairoMakie, LaTeXStrings
using SpectralMethodsTrefethen, ForwardDiff
"p11 - Chebyshev differentation of a smooth function"
function p11()
    f(x) = exp(x) * sin(5x)
    df_dx(x) = ForwardDiff.derivative(f, x)
    fig = Figure()
    for (i, N) in enumerate([10, 20])
        D, x = cheb(N)
        v = f.(x)

        Axis(fig[i, 1]; title=latexstring("\$u(x),\\,  N=$N\$"))
        scatter!(fig[i, 1], x, v)
        lines!(fig[i, 1], -1..1, f)

        error = D * v - df_dx.(x)
        Axis(fig[i, 2]; title=latexstring("error in \$u'(x),\\,  N=$N\$"))
        scatterlines!(fig[i, 2], x, error; linestyle=:dot)
    end
    return fig
end