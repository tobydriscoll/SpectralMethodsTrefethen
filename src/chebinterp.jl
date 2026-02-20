"""
    chebinterp(v)
Create a callable interpolant for the vector `v` of values given at Chebyshev nodes in [-1, 1].
"""
function chebinterp(u)
    n = length(u) - 1
    wc = [0.5; (-1).^(1:n-1); 0.5 * (-1)^n]
    xc = [cospi(k / n) for k in 0:n]
    baryinterp(xc, u, wc)
end
