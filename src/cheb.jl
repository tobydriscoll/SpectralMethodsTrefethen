"""
    cheb(N)
Chebyshev differentiation matrix and grid.
"""
function cheb(N)
    (N == 0) && return zeros(1, 1), [1.0]
    x = [sinpi(k / 2N) for k in N:-2:-N]
    c(i) = (i == 0 || i == N) ? 2 : 1
    neg1(i) = iseven(i) ? 1 : -1
    # compute the off-diagonal entries
    D = zeros(N+1, N+1)
    for i in 0:N, j in 0:N
        if i != j
            D[i+1, j+1] = neg1(i+j) * c(i) / (c(j) * (x[i+1] - x[j+1]))
        end
    end
    # negative sum trick to get the diagonal entries
    for i in 0:N
        D[i+1, i+1] = -sum(D[i+1, :])
    end
    return D, x
end