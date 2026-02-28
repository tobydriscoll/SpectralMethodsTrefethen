using CairoMakie, Printf
using LinearAlgebra, SpectralMethodsTrefethen, Polynomials
"p38 - solve u⁽⁴⁾ = exp(x), u(-1) = u(1) = u′(-1) = u′(1) = 0 (compare p13)"
function p38(N=15)
    # Construct discrete biharmonic operator:
    D, x = cheb(N)
    s = [1 / (1 - x^2) for x in x[2:N]]
    S = Diagonal([0; s; 0])
    Δ² = (Diagonal(@. 1 - x ^ 2) * D^4 - Diagonal(8x) * D^3 - 12D^2) * S
    Δ² = Δ²[2:N, 2:N]

    # Solve boundary-value problem and plot result:
    rhs = @. exp(x[2:N])
    v = Δ² \ rhs
    v = [0; v; 0]
    fig, ax, _ = scatter(x, v)
    q = chebinterp(S * v)
    u(x) = (1 - x^2) * q(x)
    lines!(ax, -1..1, u)

    # Determine exact solution and print maximum error:
    c = [1 -1 1 -1; 0 1 -2 3; 1 1 1 1; 0 1 2 3] \ exp.([-1, -1, 1, 1])
    t = Polynomial(c)
    exact(x) = exp(x) - t(x)
    maxerr = maximum(abs(u(x) - exact(x)) for x in range(-1, 1, 801))
    ax.title = @sprintf("max err = %.5g", maxerr)
    return fig
end