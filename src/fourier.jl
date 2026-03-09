"""
    fourier(N)
Fourier grid and  differentiation matrices.
"""
function fourier(N)
    @assert iseven(N) "N must be even"
    h = 2π / N
    θ = h * (1:N)
    col = [0.5 * (-1)^j * cot(j * h / 2) for j in 1:N-1]
    D = Toeplitz([0; col], [0; reverse(col)])
    c0 = -π^2 / 3h^2 - 1 / 6
    col = [c0; [ (-1)^(j+1) / 2sin(h * j / 2)^2 for j in 1:N-1 ]]
    D² = Toeplitz(col, col)
    return θ, D, D²
end
