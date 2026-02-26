using CairoMakie, LaTeXStrings, Printf, FFTW
"""
p27anim - Solve KdV eq. u_t + u u_x + u_xxx = 0 on [-π, π] by
        FFT with integrating factor v = exp(-ik³t) * û.
"""
function p27anim(N = 256, tmax = 0.006)
    @assert iseven(N) "Must choose an even value of N"
    # Set up grid and two-soliton initial data:
    x = (2π / N) * (-N/2:N/2-1)
    A, B = 25, 16
    u = @. 3A^2 * sech(0.5 * (A * (x + 2)))^2 + 3B^2 * sech(0.5 * (B * (x + 1)))^2
    û = rfft(u)

    # Set up time stepping:
    dt = 0.4 / N^2
    tplot = tmax / 150
    plotgap = ceil(Int, tplot / dt)
    dt = tplot / plotgap
    ntime = round(Int, tmax / dt)

    # Find integrating factor and solve PDE:
    k = [0:N/2-1; 0]
    ik3 = 1im * k .^ 3
    g = -0.5im * dt * k
    E = exp.(dt * ik3 / 2)
    E² = E .^ 2
    time = Observable(0.0)
    title = @lift @sprintf("t = %0.2e", $time)
    u = Observable(u)
    fig = lines(x, u;
        axis=(; xlabel=L"x", xticks=MultiplesTicks(5, π, "π"), title))
    anim = record(fig, "p27anim.mp4"; framerate=30) do io
        recordframe!(io)
        for n in 1:ntime
            a = g .* rfft(real(irfft(û, N)) .^ 2)
            b = g .* rfft(real(irfft(E .* (û + a / 2), N)) .^ 2)    # 4th-order...
            c = g .* rfft(real(irfft(E .* û + b / 2, N)) .^ 2)      # ...Runge–Kutta
            d = g .* rfft(real(irfft(E² .* û + E .* c, N)) .^ 2)
            û = @. E² * û + (E² * a + 2 * E * (b + c) + d) / 6
            if iszero(mod(n, plotgap))
                u[] = real(irfft(û, N))
                time[] = n * dt
                recordframe!(io)
            end
        end
    end
    return anim
end

p27anim()
