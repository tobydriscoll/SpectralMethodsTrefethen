using CairoMakie, LaTeXStrings
using LinearAlgebra, SpectralMethodsTrefethen
"p26 - eigenvalues of 2nd-order Chebyshev diff. matrix"
function p26(N=60)
    D, x = cheb(N)
    D² = (D^2)[2:N, 2:N]
    λ, V = eigen(-D²)

    fig = Figure(size=(600, 800))
    # Plot eigenvalues:
    M = round(maximum(λ) / N^4, sigdigits=5)
    str = "\$N = $N\$,     max \$\\!|\\lambda| = $M N^4\$"
    ax1 = Axis(fig[1, 1]; title=latexstring(str),
        xscale=log2, xlabel="mode", yscale=log10, ylabel="eigenvalue")
    scatter!(ax1, λ)
    vlines!(ax1, [2N / π], linestyle=:dash, color=:red)
    text!(ax1, 1.9N / π, 24; text=L"2\pi / N", align=(:right, :baseline))

    # Plot eigenmodes N/4 (physical) and N (nonphysical):
    vN4 = [0; V[:, Int(N / 4 - 1)]; 0]
    uN4 = chebinterp(vN4)
    ax2 = Axis(fig[2, 1]; xlabel=L"x", ylabel=L"u(x)", title=latexstring("eigenmode \$N/4\$"))
    scatter!(ax2, x, vN4)
    lines!(ax2, -1..1, uN4)

    vN = V[:, N-1];
    ax3 = Axis(fig[3, 1]; xlabel=L"x", yscale=log10, ylabel=L"|v_n|",
        title=latexstring("absolute value of eigenmode \$N\$ (log scale)"))
    scatter!(ax3, x[2:N], abs.(vN))

    rowsize!(fig.layout, 1, Relative(0.5))
    return fig
end