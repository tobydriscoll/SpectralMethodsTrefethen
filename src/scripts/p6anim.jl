using CairoMakie, LaTeXStrings, FFTW
"p6anim - variable coefficient wave equation"
function p6anim(N=128, tmax=8, Δt=π/2N)
    h = 2π / N
    x = h * (1:N)
    function fftderiv(v)
        N = length(v)
        v̂ = rfft(v)
        ŵ = im * [0:N/2-1; 0] .* v̂
        return irfft(ŵ, N)
    end

    c = @. 0.2 + sin(x - 1)^2

    # Time-stepping by leap frog formula:
    ntime = round(Int, tmax / Δt)
    Δt = tmax / ntime
    time = Observable(0.0)
    v = Observable(@. exp(-100 * (x - 1) .^ 2))
    title = @lift latexstring(@sprintf("\$t = %0.2f\$", $time))
    vold = @. exp(-100 * (x - 0.2Δt - 1)^2)
    fig = lines(x, v;
        axis=(; xlabel=L"x", xticks=MultiplesTicks(5, π, "π"), title))
    anim = record(fig, "p6anim-$N-$tmax.mp4"; framerate=60) do io
        recordframe!(io)
        for n in 1:ntime
            w = fftderiv(v[])
            vnew = vold - 2Δt * c .* w
            vold = v[]
            time[] = n * Δt
            v[] = vnew
            recordframe!(io)
        end
    end
    return anim
end