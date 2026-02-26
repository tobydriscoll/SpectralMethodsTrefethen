using CairoMakie, LaTeXStrings, SpectralMethodsTrefethen, ForwardDiff
"p18 - Chebyshev differentiation via FFT (compare p11)"
function p18()
    f(x) = exp(x) * sin(5x)
    fprime(x) = ForwardDiff.derivative(f, x)
    fig = Figure()
    for (i, N) in enumerate([10, 20])
        _, x = cheb(N)
        v = f.(x)

        Axis(fig[i, 1]; title=latexstring("\$u(x),\\,  N=$N\$"))
        scatter!(fig[i, 1], x, v)
        lines!(fig[i, 1], -1..1, f)

        error = chebfft(v) - fprime.(x)
        Axis(fig[i, 2]; title=latexstring("error in \$u'(x),\\,  N=$N\$"))
        scatter!(fig[i, 2], x, error)
        lines!(fig[i, 2], -1..1, chebinterp(error))
    end
    return fig
end