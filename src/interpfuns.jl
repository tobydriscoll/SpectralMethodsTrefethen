"""
	baryinterp(x, u)
	baryinterp(x, u, w)
Create a polynomial interpolant by the barycentric formula for the function values in
vector `u` at the node locations in vector `x`. If given, `w` is a vector of the
barycentric weights; otherwise it is computed from the nodes. The return value is
a callable function of the interpolation variable.
"""
function baryinterp(x, u, w)
    n = length(u)
    t = zeros(n)
    return function(s)
        hit = nothing
        for i in 1:n
            t[i] = w[i] / (s - x[i])
            if isinf(t[i])
                hit = i
                break
            end
        end
        isnothing(hit) ? dot(t, u) / sum(t) : u[hit]
    end
end

function baryinterp(x, u)
    n = length(u)
    w = [1 / prod(2(x[j] - x[k]) for k in 1:n if k !== j) for j in 1:n]
    baryinterp(x, u, w)
end

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

function interp2d(V, method1, method2)
    # Precompute the x interpolators.
    ux = [method1(v) for v in eachcol(V)]
    return function(x, y)
        w = [ux(x) for ux in ux]
        return method2(w)(y)
    end
end

function interp2d(V, method1, method2, x, y)
    # Precompute the x interpolators.
    Vx = reduce(hcat, [method1(v).(x) for v in eachcol(V)])
    return reduce(vcat, [method2(vx).(y') for vx in eachrow(Vx)])
end
