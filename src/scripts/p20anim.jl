using CairoMakie, LaTeXStrings, Printf
using FFTW, SpectralMethodsTrefethen
"p20anim - 2nd-order wave eq. in 2D via FFT (compare p19)"
function p20anim(N=24, tmax=1, Δt=6/N^2)
    # Grid and initial data:
    x = y = cheb(N)[2]
    DCT(v) = FFTW.r2r(v, FFTW.REDFT00)
    DST(v) = FFTW.r2r(v, FFTW.RODFT00)
    function fftderiv2(v)
        N = length(v) - 1
        x = [cospi(k / N) for k in 0:N]
        v̂ = DCT(v)
        dQ_dθ = DST(-(1:N-1) .* v̂[2:N]) / 2N    # omits n=0, n=N
        d²Q_dθ² = DCT(-(0:N).^2 .* v̂) / 2N
        w₂ = zero(v)    # return zero at boundaries
        # Chain rule for θ -> x:
        for n in 2:N
            s² = 1 / (1 - x[n]^2)
            s = sqrt(s²)
            w₂[n] = s² * (d²Q_dθ²[n] - x[n] * dQ_dθ[n-1] * s)
        end
        return w₂
    end

    V = [exp(-40 * ((x - 0.4)^2 + y^2)) for x in x, y in y]
    Vold = copy(V)
    Vxx = similar(V)
    Vyy = similar(V)

    tplot = 1 / 60
    plotgap = round(Int, tplot / Δt)
    Δt = tplot / plotgap
    ntime = round(Int, tmax / Δt)

    # Time-stepping by leap frog formula:
    time = Observable(0.0)
    title = @lift latexstring(@sprintf("\$t = %0.2f\$", $time))
    V = Observable(V)
    xx = yy = range(-1, 1, 151)
    VV = @lift interp2dgrid($V, chebinterp, chebinterp, xx, yy)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=L"x", ylabel=L"y", aspect=DataAspect(),
                limits = (-1, 1, -1, 1), title)
    heatmap!(ax, xx, yy, VV; colorrange=(-0.75, 0.75), colormap=:redsblues)
    anim = record(fig, "p20anim.mp4"; framerate=30) do io
        recordframe!(io)
        for n in 1:ntime
            Vxx = mapslices(fftderiv2, V[]; dims=[1])
            Vyy = mapslices(fftderiv2, V[]; dims=[2])
            Vnew = 2V[] - Vold + Δt^2 * (Vxx + Vyy)
            time[] = n * Δt
            Vold, V[] = V[], Vnew
            iszero(mod(n, plotgap)) && recordframe!(io)
        end
    end
    return anim
end