using CairoMakie, LaTeXStrings
using ToeplitzMatrices, SpectralMethodsTrefethen
"""
p29 - solve Poisson equation on the unit disk
         (compare p16 and p28)
"""
function p29(N = 31, M = 40)
    # Laplacian in polar coordinates:
    @assert isodd(N) && iseven(M) "Must choose odd N and even M"
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
    u = Δ \ vec(F)

    # Reshape results onto 2D grid and plot them:
    U = reshape(u, M, N2)
    U = U[[end; 1:end], :]    # repeat the periodic value for plotting
    U = [zeros(M+1) U]        # boundary values at r = 1
    X = [r * cos(θ) for θ in [0; θ], r in r[1:N2+1]]
    Y = [r * sin(θ) for θ in [0; θ], r in r[1:N2+1]]
    return contourf(X, Y, U; levels=20, colormap=:viridis,
        axis=(; xlabel=L"x", ylabel=L"y", aspect=1))
end

p29(51, 80)
