using CairoMakie, Printf, SpectralMethodsTrefethen
"""
p34anim - Allen–Cahn eq. u_t = ϵ u_xx + u - u³, u(-1) = -1, u(1) = 1\n
         (compare p6 and p32)
"""
function p34anim(N = 20, ϵ = 0.01, tmax = 1 / ϵ)
    # Differentiation matrix and initial data:
    D, x = cheb(N)
    D² = D^2                 # use full-size matrix
    D²[[1, N+1], :] .= 0     # for convenience
    v = @. 0.53x + 0.47sinpi(-1.5x)

    # Solve PDE by Euler formula and plot results:
    dt = min(0.01, 20 / (N^4 * ϵ))
    tplot = 0.2
    plotgap = ceil(Int, tplot / dt)
    dt = tplot / plotgap
    ntime = round(Int, tmax / dt)

    time = Observable(0.0)
    v = Observable(v)
    u = @lift chebinterp($v)
    str = @lift @sprintf("t = %0.1f", $time)
    fig, _ = lines(-1..1, u,
        axis=(; xlabel=L"x", ylabel=L"u(x,t)", title=str) )
    anim = record(fig, "p34anim.mp4"; framerate=60) do io
        recordframe!(io)
        for n in 1:ntime
            vnew = v[]
            vnew += dt * (ϵ * D² * (vnew - x) + vnew - vnew .^ 3)    # Euler
            v[] = vnew
            time[] = n * dt
            iszero(mod(n, plotgap)) &&recordframe!(io)
        end
    end
    return anim
end

p34anim(32)
