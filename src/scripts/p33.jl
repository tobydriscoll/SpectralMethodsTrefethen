using CairoMakie, Printf, SpectralMethodsTrefethen
"p33 - solve linear BVP u_xx = exp(4x), u′(-1) = u(1) = 0"
function p33(N = 16)
    D, x = cheb(N)
    A = D^2
    A[N+1, :] = D[N+1, :]            # Neumann condition at x = -1
    A[1, :] = I(N+1)[1, :]           # Dirichlet condition at x = 1
    f = [exp(4x) for x in x[2:N]]
    v = A \ [0; f; 0]
    fig, ax, _ = scatter(x, v)
    u = chebinterp(v)
    lines!(ax, -1..1, u)
    exact(x) = (exp(4x) - 4 * exp(-4) * (x - 1) - exp(4)) / 16
    maxerr = maximum(abs(u(x) - exact(x)) for x in range(-1, 1, 801))
    ax.title = @sprintf("max err = %0.4e", maxerr)
    return fig
end

p33()
