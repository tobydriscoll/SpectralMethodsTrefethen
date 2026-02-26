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
    r = z .- 1
    for f in [z -> 1, z -> (3 - 1 / z) / 2, z -> (23 - 16 / z + 5 / z^2) / 12]
        lines!(ax[1, 1], reim(@. r / f(z))...)
    end
    ax[1, 1].title = "Adams–Bashforth"

    # Adams–Moulton, orders 3-6:
    sfun = [z -> (5 * z + 8 - 1 / z) / 12,
        z -> (9 * z + 19 - 5 / z + 1 / z^2) / 24,
        z -> (251 * z + 646 - 264 / z + 106 / z^2 - 19 / z^3) / 720]
    for f in sfun
        lines!(ax[1, 2], reim(@. r / f(z))...)
    end
    d = @. 1 - 1 / z
    p = Polynomial([1, -1 / 2, -1 / 12, -1 / 24, -19 / 720, -3 / 160])
    lines!(ax[1, 2], reim(@. d / p(d))...)
    ax[1, 2].title = "Adams–Moulton"

    # Backward differentiation, orders 1-6:
    r = zeros(size(d))
    for i in 1:6
        r += (d .^ i) / i
        lines!(ax[2, 1], reim(r)...)
    end
    ax[2, 1].title = "backward differentiation"

    # Runge-Kutta, orders 1-4:
    taylor = [Polynomial(1 ./ factorial.(0:n)) for n in 1:4]
    w = complex(zeros(4))
    W = repeat(w, 1, length(z))
    for i in 2:length(z)
        w[1] -= (taylor[1](w[1]) - z[i])
        w[2] -= (taylor[2](w[2]) - z[i]^2) / taylor[1](w[2])
        w[3] -= (taylor[3](w[3]) - z[i]^3) / taylor[2](w[3])
        w[4] -= (taylor[4](w[4]) - z[i]^4) / taylor[3](w[4])
        W[:, i] .= w
    end
    for w in eachrow(W)
        lines!(ax[2, 2], reim(w)...)
    end
    ax[2, 2].title = "Runge–Kutta"
    return fig
end

p25()
