using CairoMakie, Printf, LaTeXStrings
using LinearAlgebra, SpectralMethodsTrefethen
"""
# p23 - eigenvalues of perturbed Laplacian on [-1,1] x [-1,1] (compare p16)
"""
function p23(N=16; perturb=false)
    # Set up tensor product Laplacian and compute 4 eigenmodes:
    D, x = cheb(N)
    y = x

    D² = (D^2)[2:N, 2:N]
    Δ = -kron(I(N-1), D²) - kron(D², I(N-1))    # Laplacian
    if perturb
        F = [exp(20 * (y - x - 1)) for x in x[2:N], y in y[2:N]]
        Δ += Diagonal(vec(F))
    end
    λ, V = eigen(Δ)

    # Reshape them to 2D grid, interpolate to finer grid, and plot:
    xx = yy = range(-1, 1, 101)
    U = zeros(N+1, N+1)
    fig = Figure(size=(600, 600))
    modes = [1 2; 3 4]
    ax = [Axis(fig[i, j]; aspect=1) for i in 1:2, j in 1:2]
    for (i, mode) in pairs(modes)
        U[2:N, 2:N] = reshape(V[:, mode], N-1, N-1)
        UU = interp2dgrid(U, chebinterp, chebinterp, xx, yy)
        UU /= maximum(abs, UU)
        str = @sprintf("λ = %0.12f π^2/4", 4λ[mode] / π^2)
        contourf!(ax[i], xx, yy, UU; levels=-1.1:0.2:1.1, colormap=:redsblues)
        ax[i].title = str
        hidedecorations!(ax[i], ticks=false)
    end
    return fig
end