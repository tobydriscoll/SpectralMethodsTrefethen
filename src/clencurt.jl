"""
    clencurt(N)
Nodes and weights for Clenshaw-Curtis quadrature.
"""
function clencurt(N)
    q = [i / N for i in 0:N]
    x = cospi.(q)
    w = zeros(N+1)
    v = ones(N-1)
    if iseven(N)
        w[1] = w[N+1] = 1 / (N^2 - 1)
        for k in 1:div(N, 2)-1
            @. v -= 2cospi(2k * q[2:N]) / (4k^2 - 1)
        end
        @. v -= cospi(N * q[2:N]) / (N^2 - 1)
    else
        w[1] = w[N+1] = 1 / N^2
        for k in 1:div(N-1, 2)
            @. v -= 2cospi(2k * q[2:N]) / (4k^2 - 1)
        end
    end
    w[2:N] .= 2v / N
    return x, w
end
