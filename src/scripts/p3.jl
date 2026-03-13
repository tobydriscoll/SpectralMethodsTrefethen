using CairoMakie
"p3 - band-limited interpolation"
function p3(h=1)
    xmax = 10
    x = -xmax:h:xmax                     # computational grid
    sinc(t) = iszero(t) ? 1.0 : sinpi(t / h) / (π * t / h)
    funs = [
        x -> float(x==0),
        x -> float(abs(x) <= 3),
        x -> max(0, 1 - abs(x)/3)
    ]
    fig = Figure()
    for (i, f) in enumerate(funs)
        v = f.(x)
        p(t) = sum( v[j] * sinc(t - x[j]) for j in eachindex(x) )
        ax = Axis(fig[i, 1], yticks=[0, 1], limits=(-xmax, xmax, -0.5, 1.5))
        scatter!(ax, x, v)
        lines!(ax, -xmax..xmax, p)
    end
    return fig
end