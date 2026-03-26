using CairoMakie, Printf, LaTeXStrings
using SpectralMethodsTrefethen, LinearAlgebra
"p37 - 2D \"wave tank\" with Neumann BCs for |y| = 1"
function p37anim(Nx=50, Ny=15, tmax=4)
    # x variable in [-A, A], Fourier:
    A = 3
    θ, _, D²θ = fourier(Nx)
    x = θ -> A * θ / π - A
    D²x = (π / A)^2 * D²θ

    # y variable in [-1, 1], Chebyshev:
    Dy, y = cheb(Ny)
    D²y = Dy^2
    idxB, idxI = [1, Ny+1], 2:Ny
    BC = -Dy[idxB, idxB] \ Dy[idxB, idxI]

    # Initial data:
    u₀(x, y) = exp(-8 * ((x + 1.5)^2 + y^2))
    V = [u₀(x, y) for x in x.(θ), y in y]
    Δt = 1 / (Nx + Ny^2)
    Vold = [u₀(x + Δt, y) for x in x.(θ), y in y]

    # Time-stepping by leap frog formula:
    ntime = ceil(Int, tmax / Δt)
    Δt = tmax / ntime
    time = Observable(0.0)
    V = Observable(V)
    θθ = range(0, 2π, 151)
    yy = range(-1, 1, 51)
    U = @lift interp2dgrid($V, fourinterp, chebinterp, θθ, yy)
    title = @lift latexstring(@sprintf("t = %0.2f", $time))
    fig = Figure()
    ax = Axis(fig[1, 1]; title, aspect=DataAspect(), xlabel=L"x", ylabel=L"y")
    heatmap!(ax, x.(θθ), yy, U;
        colormap=:redsblues, colorrange=(-0.75, 0.75), interpolate=true)
    anim = record(fig, "p37anim.mp4"; framerate=60) do io
        recordframe!(io)
        for n in 1:ntime
            Vnew = 2V[] - Vold + Δt^2 * (D²x * V[] + V[] * D²y')
            Vnew[:, idxB] .= Vnew[:, idxI] * BC'       # Neumann BCs for |y|=1
            time[] = n * Δt
            Vold, V[] = V[], Vnew
            recordframe!(io)
        end
    end
    return anim
end