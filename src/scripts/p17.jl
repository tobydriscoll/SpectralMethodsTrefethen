using CairoMakie, Printf, LaTeXStrings
using LinearAlgebra, SpectralMethodsTrefethen
"""
p17 - Helmholtz eq. u_xx + u_yy + k²u = f on [-1,1] x [-1,1] (compare p16)
"""
function p17(N=24, k=9)
    # Set up spectral grid and tensor product Helmholtz operator:
    D,x = cheb(N)
    y = x
    F = [ exp(-10*((y - 1)^2 + (x - 0.5)^2)) for x in x[2:N], y in y[2:N] ]
    D² = (D^2)[2:N,2:N]
    L = kron(I(N-1), D²) + kron(D², I(N-1)) + k^2 * I    # Helmholtz

    # Solve for u, reshape to 2D grid, and plot:
    v = L \ vec(F)
    V = zeros(N+1, N+1)
    V[2:N, 2:N] = reshape(v, N-1, N-1)
    xx = yy = range(-1, 1, 201)
    UU = interp2dgrid(V, chebinterp, chebinterp, xx, yy)

    fig = Figure(size=(600, 800))
    ax1 = Axis3(fig[1, 1]; title="Solution of Helmholtz equation",
      xlabel=L"x", ylabel=L"y", zlabel=L"u(x,y)")
    M = maximum(abs, UU)
    surface!(ax1, xx, yy, UU; colormap=:cork, colorrange=(-M, M))
    icen = div(N, 2) + 1
    text!(ax1, -1, 1, M; text="u(0, 0) = $(@sprintf("%0.11f", V[icen, icen]))")
    ax2 = Axis(fig[2, 1]; xlabel=L"x", ylabel=L"y", aspect=1)
    contour!(ax2, xx, yy, UU; levels=range(-M, M, 11), colormap=:cork)
    lines!(ax2, [-1, -1, 1, 1, -1], [-1, 1, 1, -1, -1], color=:black, linewidth=2)
    hidedecorations!(ax2; label=false, ticks=false)
    return fig
end