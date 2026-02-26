using CairoMakie, Printf, SpectralMethodsTrefethen
"p35 - Allen-Cahn eq. as in p34, but with boundary condition
    imposed explicitly (\"method (II)\")"
function p35(N = 20)
    D, x = cheb(N)
    D² = D^2      # use full-size matrix
    ϵ = 0.01
    t, Δt = 0, min(0.01, 50 / (N^4 * ϵ))
    v = @. 0.53x + 0.47sinpi(-1.5x)

    # Solve PDE by Euler formula and plot results:
    tmax = 100
    tplot = 1
    nplots = round(Int, tmax / tplot)
    plotgap = round(Int, tplot / Δt)
    Δt = tplot / plotgap
    tplots = range(0, tmax, nplots+1)
    xx = -1:0.025:1
    vv = chebinterp(v, -1, 1).(xx);
    data = [vv zeros(length(vv), nplots)]
    for i in 1:nplots
        for n in 1:plotgap
            t += Δt
            v += Δt * (ϵ * D² * (v - x) + v - v .^ 3)    # Euler
            v[1]   = 1 + sin(t / 5)^2     # BC at x=1
            v[end] = -1                   # BC at x=-1
        end
        data[:, i+1] = chebinterp(v).(xx)
    end
    fig = Figure()
    ax = Axis3(fig[1, 1]; xlabel=L"x", ylabel=L"t", zlabel=L"u(x,t)",
        limits=(nothing, nothing, (-1.05, 2)) )
    surface!(ax, xx, tplots, data;
        colormap=:redsblues, colorrange=(-1, 1))
    return fig
end