using CairoMakie, Printf, LaTeXStrings
using LinearAlgebra, SpectralMethodsTrefethen
"p16 - Poisson eq. on [-1,1] x [-1,1] with u = 0 on boundary"
function p16(N=24)
    # Set up grids and tensor product Laplacian and solve for u:
    D, x = cheb(N)
    y = x
    f(x, y) = 10sin(8x * (y - 1))
    F = [f(x, y) for x in x[2:N], y in y[2:N]]
    D² = (D^2)[2:N, 2:N]
    Δ = kron(I(N-1), D²) + kron(D², I(N-1))    # Laplacian

    fig = Figure(size=(640, 800))
    ax1 = Axis(fig[1, 1]; title="Nonzeros in the Laplacian",
      aspect=1, yreversed=true)
    spy!(ax1, Δ)
    @time v = Δ \ vec(F)    # solve problem and watch the clock

    # Reshape long 1D results onto 2D grid (reversing to usual direction):
    V = zeros(N+1, N+1)
    V[2:N, 2:N] .= reshape(v, N-1, N-1)
    icen = div(3N, 4) + 1
    value = V[icen, icen]

    # Interpolate to finer grid and plot:
    xx = yy = range(-1, 1, 201)
    UU = interp2dgrid(V, chebinterp, chebinterp, xx, yy)
    M = maximum(abs, UU)
    ax2 = Axis3(fig[2, 1]; title="Solution of Poisson equation",
      xlabel=L"x", ylabel=L"y", zlabel=L"u(x,y)")
    surface!(ax2, xx, yy, UU; colormap=:redsblues, colorrange=(-M, M))
    str = latexstring("u(2^{-1/2},2^{-1/2}) = $(@sprintf("%0.11f", value))")
    text!(ax2, -1, 1, 1.4M; text=str, align=(:left, :top))
    rowsize!(fig.layout, 1, Relative(0.35))
    return fig
end