using CairoMakie, Printf, SpectralMethodsTrefethen
"""
p34 - Allen–Cahn eq. u_t = ϵ u_xx + u - u³, u(-1) = -1, u(1) = 1\n
         (compare p6 and p32)
"""
function p34(N = 20, ϵ = 0.01, tmax = 1 / ϵ)
    # Differentiation matrix and initial data:
    D, x = cheb(N)
    D² = D^2      # use full-size matrix
    D²[[1, N+1], :] .= 0                     # for convenience
    t, dt = 0, min(0.01, 50 / (N^4 * ϵ))
    v = @. 0.53x + 0.47sinpi(-1.5x)

    # Solve PDE by Euler formula and plot results:
    tplot = 2
    nplots = round(Int, tmax / tplot)
    plotgap = round(Int, tplot / dt)
    dt = tplot / plotgap
    tplots = range(0, tmax, nplots+1)
    xx = -1:0.025:1
    vv = chebinterp(v).(xx)
    data = [vv zeros(length(vv), nplots)]
    for i in 1:nplots
        for n in 1:plotgap
            t += dt
            v += dt * (ϵ * D² * (v - x) + v - v .^ 3)    # Euler
        end
        data[:, i+1] = chebinterp(v).(xx)
    end
    fig = Figure()
    ax = Axis3(fig[1, 1]; xlabel=L"x", ylabel=L"t", zlabel=L"u(x,t)",
        limits=(nothing, nothing, (-1.05, 1.05)) )
    surface!(ax, xx, tplots, data;
        colormap=:redsblues, colorrange=(-1, 1))
    return fig
end

p34()
