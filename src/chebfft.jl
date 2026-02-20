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
