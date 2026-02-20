"""
    fourinterp(v)
Create a callable interpolant for the vector `v` of values given at equispaced (Fourier) nodes in (0, 2π].
"""
function fourinterp(v)
    N = length(v)
    @assert(iseven(N), "N must be even")
    SN(x) = sin(N * x / 2) / (N * tan(x / 2))    # δ interpolant
    return function(x)
        t = [SN(x - m * 2π / N) for m in 1:N]
        hit = isnan.(t)
        any(hit) ? v[findfirst(hit)] : dot(v, t)
    end
end
