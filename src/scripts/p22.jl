using CairoMakie, LaTeXStrings, Printf
using LinearAlgebra, ToeplitzMatrices, SpecialFunctions, SpectralMethodsTrefethen
"p22 - 5th eigenvector of Airy equation u_xx = λ*x*u"
function p22()
    fig = Figure(size=(600, 600))
    N = [12 24; 36 48]
    for (i, N) in pairs(N)
        D, x = cheb(N)
        D² = (D^2)[2:N, 2:N]
        λ, V = eigen(D², diagm(x[2:N]))         # generalized ev problem
        fifth = findfirst(>(0), λ) + 4
        λ, v = λ[fifth], V[:, fifth]
        v = [0; v; 0]                           # extend to booundary
        v = v / v[div(N, 2) + 1] * airyai(0)    # normalize
        u = chebinterp(v)
        lines(fig[i[1], i[2]], -1..1, u;
            axis=(; title=@sprintf("N = %d,  λ = %0.10f", N, λ)))
    end
    return fig
end
