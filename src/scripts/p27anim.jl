using CairoMakie, LaTeXStrings, Printf, FFTW
"""
p27anim - Solve KdV eq. u_t + u u_x + u_xxx = 0 on [-π, π] by
        FFT with integrating factor v = exp(-ik³t) * û.
"""
function p27anim(N=256, tmax=0.006, Δt=0.4/N^2)
    @assert iseven(N) "Must choose an even value of N"
    # Set up grid and two-soliton initial data:
    x = (2π / N) * (-N/2:N/2-1)
    soliton(x, a) = 3a^2 * sech(0.5 * a * x)^2
    u = @. soliton(x + 2, 25) + soliton(x + 1, 16)

    # Set up time stepping:
    tplot = tmax / 150
    plotgap = ceil(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)

    # Find integrating factor and solve PDE:
    k = [0:N/2-1; 0]
    ik³ = 1im * k .^ 3
    g = -0.5im * Δt * k
    E = exp.(Δt * ik³ / 2)
    E² = E .^ 2
    time = Observable(0.0)
    title = @lift @sprintf("t = %0.2e", $time)
    û = rfft(u)
    u = Observable(u)
    fig = lines(x, u;
        axis=(; xlabel=L"x", xticks=MultiplesTicks(5, π, "π"), title))
    anim = record(fig, "p27anim-$N-$(1000tmax).mp4"; framerate=30) do io
        recordframe!(io)
        for n in 1:ntime
            a = g .* rfft(irfft(û, N) .^ 2)
            b = g .* rfft(irfft(E .* (û + a / 2), N) .^ 2)    # 4th-order...
            c = g .* rfft(irfft(E .* û + b / 2, N) .^ 2)      # ...Runge–Kutta
            d = g .* rfft(irfft(E² .* û + E .* c, N) .^ 2)
            û = @. E² * û + (E² * a + 2 * E * (b + c) + d) / 6
            if iszero(rem(n, plotgap))
                u[] = irfft(û, N)
                time[] = n * Δt
                recordframe!(io)
            end
        end
    end
    return anim
end
