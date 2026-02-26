using CairoMakie, Printf
using SpectralMethodsTrefethen, LinearAlgebra
"p15 - solve eigenvalue BVP u_xx = λ u, u(-1) = u(1) = 0"
function p15(N = 36)
    D, x = cheb(N)
    D² = (D^2)[2:N, 2:N]
    λ, V = eigen(-D²)
    fig = Figure(size=(600, 600))
    ax = [Axis(fig[i, 1], limits=(-1, 1, -0.64, 0.64)) for i in 1:6]
    for (mode, ax) in zip(5:5:30, ax)        # plot 6 eigenvectors
        v = [0; V[:, mode]; 0]
        scatter!(ax, x, v)
        u = chebinterp(v)
        lines!(ax, -1..1, u)
        eig = "eig $mode = $(-λ[mode] * 4/π^2) π²/4"
        text!(ax, 0, 0.45; text=eig, fontsize=12, align=(:center, :baseline))
        PPW = @sprintf("%.1f", 4N / (π * mode))
        text!(ax, 0.7, 0.45; text="$PPW  ppw", fontsize=12, align=(:left, :baseline))
    end
    hidespines!.(ax)
    hidedecorations!.(ax)
    return fig
end