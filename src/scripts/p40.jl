using CairoMakie, Printf, LaTeXStrings
using LinearAlgebra, SpectralMethodsTrefethen
"p40 - eigenvalues of Orr–Sommerfeld operator (compare p38)"
function p40(Re=5772)
    N = [40 60; 80 100]
    fig = Figure(size=(640, 640))
    ax = [Axis(fig[i, j]; aspect=1, limits=(-0.8, 0.2, -1, 0)) for i in 1:2, j in 1:2]
    for (i, N) in pairs(N)
        # 2nd- and 4th-order differentiation matrices:
        D, x = cheb(N)
        D² = (D^2)[2:N, 2:N]
        S = Diagonal([0; (@. 1 / (1 - x[2:N]^2)); 0])
        D⁴ = (Diagonal(@. 1 - x^2) * D^4 - Diagonal(8x) * D^3 - 12D^2) * S
        D⁴ = D⁴[2:N, 2:N]

        # Orr–Sommerfeld operators A,B and generalized eigenvalues:
        B = D² - I
        A = (D⁴ - 2D² + I) / Re - 2im * I - 1im * Diagonal(@. 1 - x[2:N]^2) * B
        λ = eigvals(A, B)
        vlines!(ax[i], [0], color=:gray)
        scatter!(ax[i], reim(λ)...; color=real(λ), colormap=:Blues, colorrange=(-0.8, 0) )
        str = @sprintf("\$N = %d,\\; \\lambda_\\mathrm{max}\$ = %.11f", N, maximum(real, λ))
        ax[i].title = latexstring(str)
    end
    return fig
end