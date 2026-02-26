using CairoMakie, LaTeXStrings, FFTW
"""
p27 - Solve KdV eq. u_t + u u_x + u_xxx = 0 on [-π, π] by
        FFT with integrating factor v = exp(-ik³t) * û.
"""
function p27(N=256, tmax=0.006, Δt=0.4/N^2)
    @assert iseven(N) "Must choose an even value of N"
    # Set up grid and two-soliton initial data:
    x = (2π / N) * (-N/2:N/2-1)
    soliton(x, a) = 3a^2 * sech(0.5 * a * x)^2
    u = @. soliton(x + 2, 25) + soliton(x + 1, 16)

    # Set up time stepping:
    tplot = tmax / 80
    plotgap = round(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)

    # Find integrating factor and solve PDE:
    k = [0:N/2-1; 0]
    ik³ = 1im * k .^ 3
    g = -0.5im * Δt * k
    E = exp.(Δt * ik³ / 2)
    E² = E .^ 2
    û = rfft(u)
    udata = u
    tdata = [0.0]
    for n = 1:ntime
        t = n * Δt
        a = g .* rfft(irfft(û, N) .^ 2)
        b = g .* rfft(irfft(E .* (û + a / 2), N) .^ 2)    # 4th-order...
        c = g .* rfft(irfft(E .* û + b / 2, N) .^ 2)      # ...Runge–Kutta
        d = g .* rfft(irfft(E² .* û + E .* c, N) .^ 2)
        û = @. E² * û + (E² * a + 2 * E * (b + c) + d) / 6
        if iszero(rem(n, plotgap))
            u = irfft(û, N)
            udata = [udata u]
            tdata = [tdata; t]
        end
    end
    return heatmap(x, tdata, udata; colormap=:Blues, interpolate=true,
        axis=(xlabel=L"x", ylabel=L"t"))
end