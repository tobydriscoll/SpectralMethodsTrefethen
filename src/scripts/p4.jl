using CairoMakie, Printf, ToeplitzMatrices
"p4 - periodic spectral differentiation"
function p4()
    # Set up grid and differentiation matrix:
    N = 24
    h = 2π / N
    x = h * (1:N)
    col = [0.5 * (-1)^j * cot(j * h / 2) for j in 1:N-1]
    D = Toeplitz([0; col], [0; reverse(col)])

    fig = Figure()
    xticks = MultiplesTicks(5, π, "π")

    # Differentiation of a hat function:
    v = @. max(0, 1 - abs(x - π) / 2)
    Axis(fig[1, 1]; title="function", xticks)
    scatterlines!(fig[1, 1], x, v)
    Axis(fig[1, 2]; title="spectral derivative", xticks)
    scatterlines!(fig[1, 2], x, D * v)

    # Differentiation of exp(sin(x)):
    v = @. exp(sin(x))
    vprime = @. cos(x) * v
    Axis(fig[2, 1]; xticks)
    scatterlines!(fig[2, 1], x, v)
    Axis(fig[2, 2]; xticks)
    scatterlines!(fig[2, 2], x, D * v)
    str = @sprintf("max error = %0.4e", maximum(abs, D * v - vprime))
    text!(fig[2, 2], 2.2, 1, text=str, fontsize=12, align=(:left, :bottom))
    return fig
end