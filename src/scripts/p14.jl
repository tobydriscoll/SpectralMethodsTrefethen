using CairoMakie, LaTeXStrings, SpectralMethodsTrefethen
"""
p14 - solve nonlinear BVP u_xx = exp(u), u(-1) = u(1) = 0 (compare p13)
"""
function p14(N=16)
    D, x = cheb(N)
    D² = (D^2)[2:N, 2:N]
    v = zeros(N-1)
    it = 0
    while true              # fixed-point iteration
        vnew = D² \ exp.(v)
        change = maximum(abs, vnew - v)
        v = vnew
        it += 1
        (change < 1e-15 || it > 99) && break
    end
    v = [0; v; 0]
    p = chebinterp(v)
    fig = Figure()
    title = "no. steps = $it,  u(0) = $(p(0))"
    ax = Axis(fig[1, 1]; title, xlabel=L"x", ylabel=L"u(x)")
    scatter!(ax, x, v)
    lines!(ax, -1..1, p)
    return fig
end