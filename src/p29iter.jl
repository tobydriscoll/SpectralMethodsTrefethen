using CairoMakie, LaTeXStrings
using FFTW, SpectralMethodsTrefethen, LinearMaps, IterativeSolvers
"""
p29 - solve Poisson equation on the unit disk
         (compare p16 and p28)
"""
function p29(N=31, M=50)
    @assert isodd(N) && iseven(M) "Must choose odd N and even M"
    N2, M2 = div(N - 1, 2), div(M, 2)
    Dr, r = cheb(N)
    dθ = 2π / M
    θ = dθ * (1:M)

    DCT(v) = FFTW.r2r(v, FFTW.REDFT00)
    DST(v) = FFTW.r2r(v, FFTW.RODFT00)
    function fftderiv12(v)
        N = length(v) - 1
        v̂ = DCT(v)
        Q′ = DST(-(1:N-1) .* v̂[2:N]) / 2N    # omits n=0, n=N
        Q′′ = DCT(-(0:N).^2 .* v̂) / 2N
        w₁, w₂ = zero(v), zero(v)    # return zero at boundaries
        # Chain rule for θ -> x:
        for n in 2:N
            s² = 1 / (1 - r[n]^2)
            s = sqrt(s²)
            w₁[n] = -r[n] * Q′[n-1] * s
            w₂[n] = s² * (Q′′[n] - r[n] * Q′[n-1] * s)
        end
        return w₁, w₂
    end

    # Laplacian in polar coordinates:
    V = zeros(M, N+1)
    ∂rV = copy(V)
    ∂rrV = copy(V)
    ∂θθV = copy(V)
    invr = Diagonal(1 ./ r[2:N2+1])
    ξ² = (0:M2).^2
    function laplacian(v)
        V[:, 2:N2+1] = reshape(v, M, N2)
        V[:, N:-1:N2+2] = V[(@. mod1(M2 + (1:M), M)), 2:N2+1]    # extend to r ∈ [-1, 0)
        for i in axes(V, 1)
            ∂rV[i, :], ∂rrV[i, :] = fftderiv12(view(V, i, :))
        end
        for j in axes(V, 2)
            ∂θθV[:, j] = irfft(-ξ² .* rfft(V[:, j]), M)
        end
        ΔV = ∂rrV[:, 2:N2+1] + ∂rV[:, 2:N2+1] * invr + ∂θθV[:, 2:N2+1] * invr.^2
        return vec(ΔV)
    end
    Δ = LinearMap(laplacian, M * N2)

    # Right-hand side and solution for u:
    F = [-r^2 * sin(θ / 2)^4 + sin(6θ) * cos(θ / 2)^2 for θ in θ, r in r[2:N2+1]]
    u, stats = gmres(Δ, vec(F), reltol=1e-4, restart=100, log=true)
    lastres = stats[:resnorm][end]
    @info "Converged in $(stats.iters) iterations with residual norm $(lastres)"

    # Reshape results onto 2D grid and plot them:
    U = reshape(u, M, N2)
    U = U[[end; 1:end], :]    # repeat the periodic value for plotting
    U = [zeros(M+1) U]        # boundary values at r = 1
    X = [r * cos(θ) for θ in [0; θ], r in r[1:N2+1]]
    Y = [r * sin(θ) for θ in [0; θ], r in r[1:N2+1]]
    V[:, 2:N2+1] = reshape(u, M, N2)
    V[:, N:-1:N2+2] = V[(@. mod1(M2 + (1:M), M)), 2:N2+1]    # extend to r ∈ [-1, 0)
    X = [r * cos(θ) for θ in θ[1:M2+1], r in r]
    Y = [r * sin(θ) for θ in θ[1:M2+1], r in r]
    # fig, ax, plt = contourf(X, Y, U; levels=20, colormap=:amp,
    fig, ax, plt = contourf(X, Y, V[1:M2+1, :]; levels=20, colormap=:amp,
        axis=(; xlabel=L"x", ylabel=L"y", aspect=1))
    Colorbar(fig[1, 2], plt)
    return fig
end