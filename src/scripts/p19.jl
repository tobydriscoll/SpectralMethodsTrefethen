using CairoMakie, LaTeXStrings, FFTW, SpectralMethodsTrefethen
"p19 - 2nd-order wave eq. on Chebyshev grid (compare p6)"
function p19(N=80, tmax=4, Δt=8/N^2)
    _, x = cheb(N)
    v = @. exp(-200x^2)
    vold = @. exp(-200 * (x - Δt)^2)

    # Time-stepping by leap frog formula:
    tplot = 0.025
    plotgap = round(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)
    data = hcat(v, zeros(N+1, ntime))
    t = Δt * (0:ntime)
    for i = 1:ntime
        w = chebfft(chebfft(v))
        w[[1, N+1]] .= 0
        vnew = 2v - vold + Δt^2 * w
        data[:, i+1] = vnew
        vold, v = v, vnew
    end

    # Plot results:
    return heatmap(x, t[1:plotgap:end], data[:, 1:plotgap:end];
        colormap=:redsblues, colorrange=(-1, 1),
        axis=(xlabel=L"x", ylabel=L"t"))
end