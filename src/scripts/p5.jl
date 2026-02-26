using CairoMakie, Printf, FFTW
"p5 - repetition of p4 via FFT (real data only)"
function p5()
    # Differentiation of a hat function:
    N = 24
    h = 2π / N
    x = h * (1:N)
    function fftderiv(v)
        N = length(v)
        v̂ = rfft(v)
        ŵ = im * [0:N/2-1; 0] .* v̂
        return irfft(ŵ, N)
    end

    fig = Figure()
    xticks = MultiplesTicks(5, π, "π")

    # Differentiation of a hat function:
    v = @. max(0, 1 - abs(x - π) / 2)
    w = fftderiv(v)
    Axis(fig[1, 1]; title="function", xticks)
    scatterlines!(fig[1, 1], x, v)
    Axis(fig[1, 2]; title="spectral derivative", xticks)
    scatterlines!(fig[1, 2], x, w)

    # Differentiation of exp(sin(x)):
    v = @. exp(sin(x))
    vprime = @. cos(x) * v
    w = fftderiv(v)
    Axis(fig[2, 1]; xticks)
    scatterlines!(fig[2, 1], x, v)
    Axis(fig[2, 2]; xticks)
    scatterlines!(fig[2, 2], x, w)
    str = @sprintf("max error = %0.4e", maximum(abs, w - vprime))
    text!(fig[2, 2], 2.2, 1, text=str, fontsize=12, align=(:left, :bottom))
    return fig
end