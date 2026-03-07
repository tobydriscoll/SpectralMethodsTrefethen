using CairoMakie, Printf
using ToeplitzMatrices, LinearAlgebra, SpectralMethodsTrefethen
"p28a, p28b - eigenmodes of Laplacian on the disk (compare p22)"
p28() = p28a()

function p28modes(N=25, M=20)
    @assert isodd(N) && iseven(M) "Must choose odd N and even M"
    # r coordinate, ranging from -1 to 1 (N must be odd):
    N2 = div(N-1, 2)
    Dr, r = cheb(N)
    D²r = Dr^2
    D₁ = D²r[2:N2+1, 2:N2+1]
    D₂ = D²r[2:N2+1, N:-1:N2+2]
    E₁ = Dr[2:N2+1, 2:N2+1]
    E₂ = Dr[2:N2+1, N:-1:N2+2]

    # θ coordinate, ranging from 0 to 2π (M must be even):
    θ, _, D²θ = fourier(M)
    M2 = div(M, 2)

    # Laplacian in polar coordinates:
    R = Diagonal(1 ./ r[2:N2+1])
    J = [zeros(M2, M2) I; I zeros(M2, M2)]
    Δ = kron(D₁ + R * E₁, I(M)) + kron(D₂ + R * E₂, J) + kron(R^2, D²θ)

    # Compute eigenmodes:
    λ, W = eigen(-Δ, sortby=abs)
    λ = sqrt.(real(λ / λ[1]))

    # Function to generate interpolated modes:
    return function(k, θ, r)
        V = reshape(real(W[:, k]), M, N2)
        V = [zeros(M) V / maximum(abs, V)]                # BC and normalize
        V = [V V[(@. mod1(M2 + (1:M), M)), N2+1:-1:1]]    # extend to r < 0
        return λ[k], interp2dgrid(V, fourinterp, chebinterp, θ, r)
    end
end

function p28a(N=25, M=20)
    getmode = p28modes(N, M)
    # Plot eigenmodes:
    r, θ = range(-1, 1, 71), range(0, π, 111)
    X = [r * cos(θ) for θ in θ, r in r]
    Y = [r * sin(θ) for θ in θ, r in r]
    zcirc = cispi.((-200:200) / 200)
    fig = Figure(size=(640, 640))
    for (i, mode) in pairs([1 3; 6 10])
        λ, U = getmode(mode, θ, r)
        str = @sprintf("Mode %d: λ = %0.11f", mode, λ)
        ax = Axis3(fig[i[1], i[2]]; title=str, aspect=:data)
        surface!(ax, X, Y, 1.5 .+ U;
            colormap=:redsblues, colorrange=(0.5, 2.5))
        contour!(ax, X, Y, U; levels=[0], color=:darkblue, linewidth=1)
        lines!(ax, reim(zcirc)...; color=:black, linewidth=2)
        set_lights!(ax, [])
        set_ambient_light!(ax, RGBf(0.85, 0.85, 0.85))
        hidedecorations!(ax)
        hidespines!(ax)
    end
    return fig
end

function p28b(N=25, M=20)
    getmode = p28modes(N, M)
    # Plot eigenmodes:
    r, θ = range(-1, 1, 71), range(0, π, 111)
    X = [r * cos(θ) for θ in θ, r in r]
    Y = [r * sin(θ) for θ in θ, r in r]
    zcirc = cispi.((-200:200) / 200)
    fig = Figure(size=(640, 640))
    modes = reshape(1:25, 5, 5)'
    for (i, mode) in pairs(modes)
        λ, U = getmode(mode, θ, r)
        str = @sprintf("%0.4f", λ)
        ax = Axis(fig[i[1], i[2]]; title=str, titlefont=:regular, aspect=DataAspect())
        contour!(ax, X, Y, U; levels=[0], color=:darkblue, linewidth=1)
        lines!(ax, reim(zcirc)...; color=:black, linewidth=2)
        hidedecorations!(ax)
        hidespines!(ax)
    end
    return fig
end