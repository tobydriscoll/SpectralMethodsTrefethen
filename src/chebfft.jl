DCT(v) = FFTW.r2r(v, FFTW.REDFT00)
DST(v) = FFTW.r2r(v, FFTW.RODFT00)

"""
    chebfft(v)
Differentiate values given at Chebyshev points via the FFT.
"""
function chebfft(v)
    N = length(v) - 1
    (N == 0) && return zero(v)
    x = [cospi(k / N) for k in 0:N]
    v̂ = DCT(v)
    w = similar(v)
    if N > 1
        ξ = 1:N-1
        s = @. -sqrt(1 - x[2:N]^2)    # chain rule denominator
        w[2:N] = DST(-ξ .* v̂[2:N]) ./ (2N * s)
    end
    w[1] = N/2 * v̂[N+1] + sum(i^2 * v̂[i+1] for i in 0:N-1) / N
    w[N+1] = (-1)^(N+1) * N/2 * v̂[N+1] +
        sum((-1)^(i+1) * i^2 * v̂[i+1] for i in 0:N-1) / N
    return w
end
