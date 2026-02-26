using  CairoMakie, LinearAlgebra, SpectralMethodsTrefethen
"p24 - pseudospectra of Davies's complex harmonic oscillator"
function p24(N = 70, L = 6)
    # Eigenvalues:
    D, x = cheb(N)
    x, D = L * x, D / L            # rescale to [-L,L]
    A = -D^2 + (1 + 3im) * Diagonal(x .^ 2)
    A = A[2:N, 2:N]
    λ = eigvals(A)

    # Pseudospectra:
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=L"\mathrm{Re}(z)", ylabel=L"\mathrm{Im}(z)",
        limits=(0, 50, 0, 40), aspect=DataAspect())
    x, y = 0:0.5:50, 0:0.5:40
    minsigma(z) = minimum(svdvals(z * I - A))
    σ_min = [minsigma(complex(x, y)) for x in x, y in y]
    cplt = contourf!(x, y, log10.(σ_min); levels=-5:0.5:-0.5)
    scatter!(ax, reim(λ)...; color=:black, markersize=6)
    Colorbar(fig[1, 2], cplt)
    return fig
end

p24()
