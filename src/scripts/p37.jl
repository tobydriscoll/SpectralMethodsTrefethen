using CairoMakie, Printf, LaTeXStrings
using SpectralMethodsTrefethen, LinearAlgebra, ToeplitzMatrices
"p37 - 2D \"wave tank\" with Neumann BCs for |y| = 1"
function p37(Nx = 50, Ny = 15)
    # x variable in [-A, A], Fourier:
    A = 3
    dθ = 2π / Nx
    θ = dθ * (1:Nx)
    c0 = -π^2 / 3dθ^2 - 1 / 6
    col = [0.5 * (-1)^(k + 1) / sin(k * dθ / 2)^2 for k in 1:Nx-1]
    D²θ = Toeplitz([c0; col], [c0; col])
    x = θ -> A * θ / π - A
    D²x = (π / A)^2 * D²θ

    # y variable in [-1, 1], Chebyshev:
    Dy, y = cheb(Ny)
    D²y = Dy^2
    BC = -Dy[[1, Ny+1], [1, Ny+1]] \ Dy[[1, Ny+1], 2:Ny]

    # Initial data:
    V = [exp(-8 * ((x + 1.5)^2 + y^2)) for x in x.(θ), y in y]
    dt = 1 / (Nx + Ny^2)
    Vold = [exp(-8 * ((x + dt + 1.5)^2 + y^2)) for x in x.(θ), y in y]

    # Time-stepping by leap frog formula:
    tmax = 4
    ntime = ceil(Int, tmax / dt)
    dt = tmax / ntime
    time = Observable(0.0)
    V = Observable(V)
    θθ = range(0, 2π, 111)
    yy = range(-1, 1, 51)
    VV = @lift interp2d($V, fourinterp, chebinterp, θθ, yy)
    title = @lift latexstring(@sprintf("t = %0.2f", $time))
    fig = Figure()
    ax = Axis(fig[1, 1]; title, aspect=DataAspect(), xlabel=L"x", ylabel=L"y")
    heatmap!(ax, x.(θθ), yy, VV; colormap=:redsblues, colorrange=(-0.75, 0.75), interpolate=true)
    anim = record(fig, "p37anim.mp4"; framerate=60) do io
        recordframe!(io)
        for i in 1:ntime
            Vnew = 2V[] - Vold + dt^2 * (D²x * V[] + V[] * D²y')
            Vnew[:, [1, Ny+1]] .= Vnew[:, 2:Ny] * BC'       # Neumann BCs for |y|=1
            time[] = i * dt
            Vold, V[] = V[], Vnew
            recordframe!(io)
        end
    end
    return anim
end