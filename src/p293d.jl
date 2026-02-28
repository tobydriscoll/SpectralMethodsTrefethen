using CairoMakie, LaTeXStrings
using FFTW, SpectralMethodsTrefethen, LinearMaps, IterativeSolvers
"""
p29 - solve Poisson equation on the unit disk
         (compare p16 and p28)
"""
function p29(N=[31, 50, 20])
    # @assert isodd(N) && iseven(M) "Must choose odd N and even M"
    N2 = [div(N[1] - 1, 2), div(N[2], 2)]
    Dr, r = cheb(N[1])
    dθ = 2π / N[2]
    θ = dθ * (1:N[2])
    Dz, z = cheb(N[3])

    DCT(v) = FFTW.r2r(v, FFTW.REDFT00)
    DST(v) = FFTW.r2r(v, FFTW.RODFT00)
    function fftderiv12(v)
        M = length(v) - 1
        v̂ = DCT(v)
        Q′ = DST(-(1:M) .* v̂[2:M+1]) / 2M    # omits n=0, n=M
        Q′′ = DCT(-(0:M).^2 .* v̂) / 2M
        w₁, w₂ = zero(v), zero(v)    # return zero at boundaries
        # Chain rule for θ -> x:
        for n in 2:M+1
            s² = 1 / (1 - r[n]^2)
            s = sqrt(s²)
            w₁[n] = -r[n] * Q′[n-1] * s
            w₂[n] = s² * (Q′′[n] - r[n] * Q′[n-1] * s)
        end
        return w₁, w₂
    end

    # Laplacian in polar coordinates:
    V = zeros(N[1]+1, N[2], N[3]+1)
    idx = CartesianIndices((2:N2[1]+1, 1:N[2], 2:N[3]))
    ΔV = similar(V, size(idx))
    ∂rV = copy(V)
    ∂rrV = copy(V)
    ∂θθV = copy(V)
    ∂zzV = copy(V)
    invr = 1 ./ r[2:N[1]]
    ξ² = (0:N2[2]).^2
    function laplacian(v)
        V[idx] = reshape(v, N2[1], N[2], N[3]-1)    # shape to the grid
        V[N[1]:-1:N2[1]+2, :, :] = V[2:N2[1]+1, (@. mod1(N2[2] + (1:N[2]), N[2])), :]     # extend to r ∈ [-1, 0)
        for (i, v) in pairs(eachslice(V, dims=(1,3)))
            ∂θθV[i[1], :, i[2]] = irfft(-ξ² .* rfft(v), N[2])
        end
        for (i, v) in pairs(eachslice(V, dims=(2,3)))
            ∂rV[:, i], ∂rrV[:, i] = fftderiv12(v)
            ∂rV[2:N[1], i] .*= invr
            ∂θθV[2:N[1], i] .*= invr.^2
        end
        for (i, v) in pairs(eachslice(V, dims=(1,2)))
            ∂zzV[i, :] = fftderiv12(v)[2]
        end
        @. ΔV = ∂rrV[idx] + ∂rV[idx] + ∂θθV[idx] + ∂zzV[idx]
        return vec(ΔV)
    end
    Δ = LinearMap(laplacian, N2[1] * N[2] * (N[3]-1))

    # Right-hand side and solution for u:
    F = [-r^2 * sin(θ / 2)^4 + sin(6θ) * cos(θ / 2)^2 + z for r in r[2:N2[1]+1], θ in θ, z in z[2:N[3]]]
    Main.@exfiltrate
    u, stats = gmres(Δ, vec(F), reltol=1e-4, maxiter=500, restart=60, log=true, verbose=true)
    lastres = stats[:resnorm][end]
    @info "Converged in $(stats.iters) iterations with residual norm $(lastres)"

    # Reshape results onto 2D grid and plot them:
    # U = reshape(u, N2[1], N[2], N[3]-1)
    # U = U[[end; 1:end], :]    # repeat the periodic value for plotting
    # U = [zeros(M+1) U]        # boundary values at r = 1
    # X = [r * cos(θ) for θ in [0; θ], r in r[1:N2+1]]
    # Y = [r * sin(θ) for θ in [0; θ], r in r[1:N2+1]]
    # V[:, 2:N2+1] = reshape(u, M, N2)
    # V[:, N:-1:N2+2] = V[(@. mod1(M2 + (1:M), M)), 2:N2+1]    # extend to r ∈ [-1, 0)
    # X = [r * cos(θ) for θ in θ[1:M2+1], r in r]
    # Y = [r * sin(θ) for θ in θ[1:M2+1], r in r]
    # # fig, ax, plt = contourf(X, Y, U; levels=20, colormap=:amp,
    # fig, ax, plt = contourf(X, Y, V[1:M2+1, :]; levels=20, colormap=:amp,
    #     axis=(; xlabel=L"x", ylabel=L"y", aspect=1))
    # Colorbar(fig[1, 2], plt)
    return u
end