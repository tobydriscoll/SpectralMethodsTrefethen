"""
    cheb(N)
Chebyshev differentiation matrix and grid.
"""
function cheb(N)
    (N == 0) && return zeros(1, 1), [1.0]
    x = [sinpi(k / N) for k in N:-2:-N]
    c = [2; ones(N-1); 2]
    @. c[2:2:end] = -c[2:2:end]
    # compute the off-diagonal entries
    D = [c[i] / (c[j] * (x[i] - x[j] + (i==j))) for i in eachindex(x), j in eachindex(x)]
    # negative sum trick to get the diagonal entries
    for i in eachindex(x)
        D[i, i] = -sum(D[i, j] for j in eachindex(x) if j != i)
    end
    return D, x
end
