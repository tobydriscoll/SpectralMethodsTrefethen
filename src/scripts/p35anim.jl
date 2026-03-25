using CairoMakie, Printf, SpectralMethodsTrefethen
"""
p35anim - Allen-Cahn eq. as in p34, but with boundary condition
              imposed explicitly ("method (II)") w/animated solution
"""
function p35anim(N=20, ϵ=0.01, tmax=1/2ϵ)
    # Differentiation matrix and initial data:
    D, x = cheb(N)
    D² = D^2     # use full-size matrix
    v = @. 0.53 * x + 0.47 * sin(-1.5π * x)

    # Solve PDE by Euler formula and plot results:
    Δt = min(0.01, 50 / (N^4 * ϵ))
    tplot = 0.2
    plotgap = ceil(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)

    time = Observable(0.0)
    v = Observable(v)
    u = @lift chebinterp($v)
    str = @lift @sprintf("t = %0.1f", $time)
    fig, ax, _ = lines(-1..1, u,
        axis=(; xlabel=L"x", ylabel=L"u(x,t)", title=str, limits=(-1, 1, -1, 2)))
    anim = record(fig, "p35anim.mp4"; framerate=15) do io
        recordframe!(io)
        for n in 1:ntime
            time[] = n * Δt
            vnew = v[]
            vnew += Δt * (ϵ * D² * (vnew - x) + vnew - vnew .^ 3)    # Euler
            vnew[1]   = 1 + sin(time[] / 5)^2    # BC at x=1
            vnew[end] = -1                       # BC at x=-1
            v[] = vnew
            iszero(mod(n, plotgap)) && recordframe!(io)
        end
    end
    return anim
end