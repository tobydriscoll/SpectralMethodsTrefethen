using CairoMakie, Printf, LinearAlgebra, SpectralMethodsTrefethen
"p13 - solve linear BVP u_xx = exp(4x), u(-1) = u(1) = 0"
function p13(N = 16)
    D, x = cheb(N)
    D² = D^2
    D² = D²[2:N, 2:N]                    # boundary conditions
    rhs = @. exp(4x[2:N])
    v = D² \ rhs                         # Poisson eq. solved here
    v = [0; v; 0]
    xx = range(-1, 1, 801)
    u = chebinterp(v)                    # interpolate grid data
    exact(x) = @. (exp(4x) - sinh(4) * x - cosh(4)) / 16;
    maxerr = maximum(abs, u.(xx) - exact.(xx))
    title = @sprintf("max err = %.3g", maxerr)
    fig = lines(-1..1, u; axis=(; title, xlabel=L"x", ylabel=L"u(x)"))
    scatter!(x, v)
    fig
end