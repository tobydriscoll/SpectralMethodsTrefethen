using CairoMakie, Printf, SpectralMethodsTrefethen
"p32 - solve u_xx = exp(4x), u(-1) = 0, u(1) = 1 (compare p13)"
function p32(N = 16)
    D, x = cheb(N)
    D² = (D^2)[2:N, 2:N]                   # boundary conditions
    f = [exp(4x) for x in x[2:N]]
    v = D² \ f                           # Poisson eq. solved here
    v = [0; v; 0] + (x .+ 1) / 2
    fig, ax, _ = scatter(x, v)
    u = chebinterp(v)
    lines!(ax, -1..1, u)
    exact(x) = (exp(4 * x) - sinh(4) * x - cosh(4)) / 16 + (x + 1) / 2
    maxerr = maximum(abs(u(x) - exact(x)) for x in range(-1, 1, 801))
    ax.title = @sprintf("max err = %0.4e", maxerr)
    return fig
end

p32()
