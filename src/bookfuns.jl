## Implmentations of functions called throughout the book.

"""
    cheb(N)
Chebyshev differentiation matrix and grid.
"""
function cheb(N)
    N == 0 && return zeros(1, 1), [1]
    x = [cospi(k / N) for k in 0:N]
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
# end of cheb()

function fftdiff(V::AbstractVector{<:Real}, N::Integer)
    V̂ = rfft(V)
    ξ = [0:N-1; 0]
    return real(V̂), irfft(im * ξ .* V̂, 2N)
end

function fftdiff(V::AbstractVector{<:Complex}, N::Integer)
    V̂ = fft(V)
    ξ = [0:N-1; 0; -N+1:-1]
    return V̂, ifft(im * ξ .* V̂)
end
# end of fftdiff()

"""
    chebfft(v)
Differentiate values given at Chebyshev points via the FFT.
"""
function chebfft(v)
    # This can be rewritten to use dct in order to exploit symmetry more efficiently.
    N = length(v) - 1
    N == 0 && return [0]
    x = [cospi(k / N) for k in 0:N]
    V = [v; reverse(v[2:N])]              # transform x -> theta
    V̂, W = fftdiff(V, N)
    w = similar(v)
    @. w[2:N] = -W[2:N] / sqrt(1 - x[2:N]^2)      # transform theta -> x
    w[1] = N/2 * V̂[N+1] + sum(i^2 * V̂[i+1] for i in 0:N-1) / N
    w[N+1] = (-1)^(N+1) * N/2 * V̂[N+1] + sum((-1)^(i+1) * i^2 * V̂[i+1] for i in 0:N-1) / N
    return w
end
# end of chebfft()

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
# end of clencurt()

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
# end of gauss()
