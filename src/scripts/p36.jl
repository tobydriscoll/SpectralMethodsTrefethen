using CairoMakie, Printf, LaTeXStrings, SpectralMethodsTrefethen, LinearAlgebra
"p36 - Laplace eq. on [-1,1] x [-1,1] with nonzero BCs"
function p36(N = 24)
    # Set up grid and 2D Laplacian, boundary points included:
    D, x = cheb(N)
    y = x
    D² = D^2
    Δ = kron(I(N+1), D²) + kron(D², I(N+1))

    # Impose boundary conditions by replacing appropriate rows of L:
    isboundary(x, y) = (abs(x) == 1 || abs(y) == 1)
    b = vec([isboundary(x, y) for x in x, y in y])       # boundary pts
    Δ[b, :] .= 0
    Δ[b, b] .= I(4N)
    rhs = zeros((N+1)^2)
    g(x, y) = (y == 1) && (x < 0) ? sinpi(x)^4 : (x == 1) ? 0.2sinpi(3y) : 0
    rhs[b] .= [g(x, y) for x in x, y in y if isboundary(x, y)]

    # Solve Laplace equation, reshape to 2D, and plot:
    v = Δ \ rhs
    V = reshape(v, N+1, N+1)
    xx = yy = range(-1, 1, 121)
    UU = interp2d(V, chebinterp, chebinterp, xx ,yy)
    u = interp2d(V, chebinterp, chebinterp)
    title = latexstring(@sprintf("u(0, 0) = %0.10f", u(0, 0)))
    fig = Figure()
    ax = Axis3(fig[1, 1]; title, xlabel=L"x", ylabel=L"y", zlabel=L"u(x,y)")
    surface!(ax, xx, yy, UU)
    return fig
end