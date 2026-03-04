using CairoMakie, Polynomials
"p25 - stability regions for ODE formulas"
function p25()
    fig = Figure(size=(640, 640))
    ax = [Axis(fig[i, j]; aspect=DataAspect()) for i in 1:2, j in 1:2]
    for ax in ax
        hlines!(ax, [0]; color=:gray)
        vlines!(ax, [0]; color=:gray)
    end

    z = [cispi(x) for x in (0:200) / 100]
    # Adams–Bashforth, orders 1-3:
    ρz = @. z - 1
    for σ in Polynomial.([[1], [-1, 3] / 2, [5, -16, 23] / 12])
        w = @. ρz / σ(z)
        lines!(ax[1, 1], reim(w)...)
        @. ρz *= z
    end
    ax[1, 1].title = "Adams–Bashforth"

    # Adams–Moulton, orders 3-6:
    ρz = @. z^2 - z
    for σ in Polynomial.([
            [-1, 8, 5] / 12,
            [1, -5, 19, 9] / 24,
            [-19, 106, -264, 646, 251] / 720,
            [27, -173, 482, -798, 1427, 475] / 1440
        ])
        w = @. ρz / σ(z)
        lines!(ax[1, 2], reim(w)...)
        @. ρz *= z
    end
    ax[1, 2].title = "Adams–Moulton"

    # Backward differentiation, orders 1-6:
    d = @. 1 - 1 / z
    r = zero(d)
    for n in 1:6
        @. r += (d ^ n) / n
        lines!(ax[2, 1], reim(r)...)
    end
    ax[2, 1].title = "backward differentiation"

    # Runge-Kutta, orders 1-4:
    T = [Polynomial(1 ./ factorial.(0:n)) for n in 0:4]
    w = zero(z)
    for n in 1:4
        for i in 1:length(z)-1
            w[i+1] = w[i] - (T[n+1](w[i]) - z[i+1]^n) / T[n](w[i])
        end
        lines!(ax[2, 2], reim(w)...)
    end
    ax[2, 2].title = "Runge–Kutta"
    return fig
end