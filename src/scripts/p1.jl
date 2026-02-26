using CairoMakie, LaTeXStrings
using ToeplitzMatrices
"p1 - convergence of fourth-order finite differences"
function p1()
    # For various N, set up grid in [-π, π] and function u(x):
    N = [2^j for j in 3:12]
    errors = []
    for N in N
        h = 2π / N
        x = [-π + j * h for j in 1:N]
        u = @. exp(sin(x)^2)
        uprime = @. 2sin(x) * cos(x) * u

        # Construct sparse fourth-order differentiation matrix:
        col = zeros(N)
        col[[2, 3, N-1, N]] = [-2/3, 1/12, -1/12, 2/3] / h
        D = Toeplitz(col, -col)

        # Compute max of abs(D*u - uprime):
        error = maximum(abs, D * u - uprime)
        push!(errors, error)
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; title="Convergence of fourth-order finite differences",
        xscale=log10, yscale=log10, xlabel=L"N", ylabel="error")
    scatter!(ax, N, errors)
    lines!(ax, N, 1 ./ N .^ 4; linestyle=:dash, color=:gray)
    text!(ax, 105, 1e-8; text=L"N^{-4}")
    return fig
end