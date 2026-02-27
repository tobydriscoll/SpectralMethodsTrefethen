using CairoMakie, LaTeXStrings
using ToeplitzMatrices, SpectralMethodsTrefethen
"p29 - solve Poisson equation on the unit disk (compare p16 and p28)"
function p29(N=31, M=40)
    @assert isodd(N) && iseven(M) "Must choose odd N and even M"
    # Set up r and θ:
    N2, M2 = div(N - 1, 2), div(M, 2)
    Dr, r = cheb(N)
    D²r = Dr^2
    D₁ = D²r[2:N2+1, 2:N2+1]
    D₂ = D²r[2:N2+1, N:-1:N2+2]
    E₁ = Dr[2:N2+1, 2:N2+1]
    E₂ = Dr[2:N2+1, N:-1:N2+2]
    dθ = 2π / M
    θ = dθ * (1:M)
    c0 = -π^2 / 3dθ^2 - 1 / 6
    col = [0.5 * (-1)^(k + 1) / sin(k * dθ / 2)^2 for k in 1:M-1]
    D²θ = Toeplitz([c0; col], [c0; col])

    # Laplacian in polar coordinates:
    R = Diagonal(1 ./ r[2:N2+1])
    J = [zeros(M2, M2) I; I zeros(M2, M2)]
    Δ = kron(D₁ + R * E₁, I(M)) + kron(D₂ + R * E₂, J) + kron(R^2, D²θ)

    # Right-hand side and solution for u:
    F = [-r^2 * sin(θ / 2)^4 + sin(6θ) * cos(θ / 2)^2 for θ in θ, r in r[2:N2+1]]
    v = Δ \ vec(F)

    # Reshape results onto 2D grid and plot them:
    V = reshape(v, M, N2)
    V = [zeros(M) V]        # boundary values at r = 1
    V = [V V[(@. mod1(M2 + (1:M), M)), N2+1:-1:1]]    # extend to r ∈ [-1, 0)
    rr, θθ = range(-1, 1, 71), range(0, π, 111)
    U = interp2dgrid(V, fourinterp, chebinterp, θθ, rr)
    X = [r * cos(θ) for θ in θθ, r in rr]
    Y = [r * sin(θ) for θ in θθ, r in rr]
    fig, ax, plt = contourf(X, Y, U; levels=20, colormap=:amp,
        axis=(; xlabel=L"x", ylabel=L"y", aspect=1))
    Colorbar(fig[1, 2], plt)
    return fig
end