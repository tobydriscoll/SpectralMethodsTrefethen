"""
    fourinterp(v)
Create a callable interpolant for the vector `v` of values given at equispaced (Fourier)
nodes in (0, 2π].
"""
function fourinterp(v)
    N = length(v)
    @assert(iseven(N), "N must be even")
    h = 2π / N
    # Periodic sinc function:
    SN(x) = iszero(mod(x, 2π)) ? 1.0 : sin(N * x / 2) / (N * tan(x / 2))
    return x -> sum(v[j] * SN(x - j * h) for j in 1:N)
end