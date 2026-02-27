"""
    gauss(N)
Nodes and weights for Gauss quadrature.
"""
function gauss(N)
    β = [0.5 / sqrt(1 - 1 / 4i^2) for i in 1:N-1]
    T = SymTridiagonal(zeros(N), β)
    x, V = eigen(T)
    return x, 2V[1, :] .^ 2
end