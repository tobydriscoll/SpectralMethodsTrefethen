using CairoMakie, LaTeXStrings, LinearAlgebra, ToeplitzMatrices
"""
p21 - eigenvalues of Mathieu operator -u_xx + 2 q cos(2x) u\n
         (compare p8 and p. 724 of Abramowitz & Stegun)
"""
function p21(N = 42)
    h = 2π / N
    x = h * (1:N)
    c0 = -π^2 / 3h^2 - 1 / 6
    col = [c0; [ (-1)^(j+1) / 2sin(h * j / 2)^2 for j in 1:N-1 ]]
    D² = Toeplitz(col, col)    # second-order differentiation matrix
    q = 0:.2:15
    data = zeros(11, length(q))
    for (i, q) in enumerate(q)
        λ = eigvals( -D² + 2q * Diagonal(cos.(2x)) )
        data[:, i] = λ[1:11]
    end
    fig = Figure(size=(320, 500))
    ax = Axis(fig[1, 1]; xlabel=L"q", ylabel=L"λ", yticks=-24:4:32,
        limits=(0, 15, -24, 32))
    series!(ax, q, data[1:2:end, :], color=:atlantic)
    series!(ax, q, data[2:2:end, :], linestyle=:dot, color=:atlantic)
    return fig
end