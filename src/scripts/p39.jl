using CairoMakie, Printf
using LinearAlgebra, SpectralMethodsTrefethen, Polynomials
"p39 - eigenmodes of biharmonic on a square with clamped BCs (compare p38)"
function p39(N = 17)
    # Construct spectral approximation to biharmonic operator:
    D, x = cheb(N)
    y = x
    D² = (D^2)[2:N, 2:N]
    s = [1 / (1 - x^2) for x in x[2:N]]
    S = Diagonal([0; s; 0])
    D⁴ = (Diagonal(@. 1 - x ^ 2) * D^4 - Diagonal(8x) * D^3 - 12D^2) * S
    D⁴ = D⁴[2:N, 2:N]
    Δ² = kron(I(N-1), D⁴) + kron(D⁴, I(N-1)) + 2kron(D², I(N-1)) * kron(I(N-1), D²)

    # Find and plot 25 eigenmodes:
    λ, V = eigen(Δ²)
    λ, V = real(λ[1:25]), real(V[:, 1:25])
    λ = sqrt.(λ / λ[1])
    xx = yy = -1:0.01:1
    square = [1 + 1im, -1 + 1im, -1 - 1im, 1 - 1im, 1 + 1im]
    fig = Figure(size=(640, 640))
    modes = reshape(1:25, 5, 5)'
    for (i, mode) in pairs(modes)
        str = @sprintf("%0.6g", λ[mode])
        U = zeros(N+1, N+1)
        U[2:N, 2:N] = reshape(V[:, mode], N-1, N-1)
        UU = interp2d(U, chebinterp, chebinterp, xx, yy)
        ax = Axis(fig[i[1], i[2]]; title=str, titlefont=:regular, aspect=DataAspect())
        contour!(ax, xx, yy, UU; levels=[0], color=:black, linewidth=1)
        lines!(ax, reim(square)...; linewidth=1, color=:black)
        hidedecorations!(ax)
        hidespines!(ax)
    end
    return fig
end