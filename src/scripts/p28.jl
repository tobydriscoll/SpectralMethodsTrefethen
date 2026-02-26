using CairoMakie, Printf
using ToeplitzMatrices, LinearAlgebra, SpectralMethodsTrefethen
"p28 - eigenmodes of Laplacian on the disk (compare p22)"
function p28(N = 25, M = 20; version=:a, interpolate=false)
    @assert isodd(N) && iseven(M) "Must choose odd N and even M"
    # r coordinate, ranging from -1 to 1 (N must be odd):
    N2 = div(N-1, 2)
    Dr, r = cheb(N)
    D²r = Dr^2
    D₁ = D²r[2:N2+1, 2:N2+1]
    D₂ = D²r[2:N2+1, N:-1:N2+2]
    E₁ = Dr[2:N2+1, 2:N2+1]
    E₂ = Dr[2:N2+1, N:-1:N2+2]

    # t = theta coordinate, ranging from 0 to 2π (M must be even):
    dθ = 2π / M
    θ = dθ * (1:M)
    M2 = div(M, 2)
    c0 = -π^2 / 3dθ^2 - 1 / 6
    col = [0.5 * (-1)^(k + 1) / sin(k * dθ / 2)^2 for k in 1:M-1]
    D²θ = Toeplitz([c0; col], [c0; col])

    # Laplacian in polar coordinates:
    R = Diagonal(1 ./ r[2:N2+1])
    J = [zeros(M2, M2) I; I zeros(M2, M2)]
    Δ = kron(D₁ + R * E₁, I(M)) + kron(D₂ + R * E₂, J) + kron(R^2, D²θ)

    # Compute four eigenmodes:
    λ, V = eigen(-Δ, sortby=abs)
    λ = sqrt.(real(λ / λ[1]))

    function getmode(V, k)
        U = reshape(real(V[:, k]), M, N2)
        U = [zeros(M) U / maximum(abs, U)]
        return U[[M; 1:M], :]
    end

    function getmode(V, k, θ, r)
        U = reshape(real(V[:, k]), M, N2)     # shape to the grid
        U = [zeros(M) U / maximum(abs, U)]    # BC and normalize
        U = [U U[(@. mod1(M2 + (1:M), M)), N2:-1:1]]    # extend to r ∈ [-1, 0)
        U = interp2d(U, fourinterp, chebinterp, θ, r)
        return U
    end

    # Plot eigenmodes:
    if interpolate
        θ = range(0, 2π, 101)
        r = range(0, 1, 51)
    else
        θ = [0; θ]
        r = r[1:N2+1]
    end
    X = [r * cos(θ) for θ in θ, r in r]
    Y = [r * sin(θ) for θ in θ, r in r]
    zcirc = cispi.((-200:200) / 200)
    fig = Figure(size=(640, 640))
    if version == :a
        modes = [1 3; 6 10]
        for (i, mode) in pairs(modes)
            str = @sprintf("Mode %d: λ = %0.11f", mode, λ[mode])
            U = interpolate ? getmode(V, mode, θ, r) : getmode(V, mode)
            ax = Axis3(fig[i[1], i[2]]; title=str, aspect=:data)
            surface!(ax, X, Y, 1.5 .+ U;
                colormap=:redsblues, colorrange=(0.5, 2.5))
            contour!(ax, X, Y, U; levels=[0], color=:black, linewidth=1.5)
            lines!(ax, reim(zcirc)...; color=:black, linewidth=1.5)
            hidedecorations!(ax)
            hidespines!(ax)
        end
    elseif version == :b
        modes = reshape(1:25, 5, 5)'
        for (i, mode) in pairs(modes)
            str = @sprintf("%0.4f", λ[mode])
            U = interpolate ? getmode(V, mode, θ, r) : getmode(V, mode)
            ax = Axis(fig[i[1], i[2]]; title=str, titlefont=:regular, aspect=DataAspect())
            contour!(ax, X, Y, U; levels=[0], color=:black, linewidth=1)
            lines!(ax, reim(zcirc)...; color=:black, linewidth=1)
            hidedecorations!(ax)
            hidespines!(ax)
        end
    end
    return fig
end

p28()
