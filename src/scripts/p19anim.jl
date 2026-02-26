using CairoMakie, Printf, LaTeXStrings
using FFTW, SpectralMethodsTrefethen
"p19anim - 2nd-order wave eq. on Chebyshev grid (compare p6)"
function p19anim(N=80, tmax=4)
    _, x = cheb(N)
    Δt = 8 / N^2

    tplot = 0.005
    plotgap = round(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)

    # Time-stepping by leap frog formula:
    time = Observable(0.0)
    title = @lift latexstring(@sprintf("\$t = %0.2f\$", $time))
    v = Observable(@. exp(-200x^2))
    vold = @. exp(-200 * (x - Δt)^2)
    fig = lines(x, v; axis=(; xlabel=L"x", title, limits=(-1, 1, -1, 1)))
    anim = record(fig, "p19anim.mp4"; framerate=60) do io
        recordframe!(io)
        for n in 1:ntime
            w = chebfft(chebfft(v[]))
            w[[1, N+1]] .= 0
            vnew = 2v[] - vold + Δt^2 * w
            vold = v[]
            v[] = vnew
            time[] = n * Δt
            iszero(mod(n, plotgap)) && recordframe!(io)
        end
      end
    return anim
end