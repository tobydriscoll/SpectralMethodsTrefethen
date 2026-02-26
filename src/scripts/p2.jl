using CairoMakie, LaTeXStrings, ToeplitzMatrices
"p2 - convergence of periodic spectral method (compare p1)"
function p2()
    # For various N (even), set up grid as before:
    N = 2:2:100
    errors = zeros(size(N))
    for (k, N) in enumerate(N)
        h = 2π / N
        x = [-π + j * h for j in 1:N]
        u = @. exp(sin(x))
        uprime = @. cos(x) * u

        # Construct spectral differentiation matrix:
        col = [0.5 * (-1)^i * cot(i * h / 2) for i in 1:N-1]
        D = Toeplitz([0; col], [0; reverse(col)])

        # Compute max(abs(D*u - uprime)):
        error = maximum(abs, D * u - uprime)
        errors[k] = error
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; title="Convergence of spectral differentiation",
        xscale=log10, yscale=log10, xlabel=L"N", ylabel="error")
    scatter!(ax, N, errors)
    return fig
end