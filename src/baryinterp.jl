"""
    baryweights(x)
Compute the barycentric weights for the nodes in vector `x`. Returns a vector of the
same length as `x`.
"""
function baryweights(x)
    n = length(x)
    a, b = extrema(x)
    x̃ = 4x / (b - a)    # rescale to avoid overflow/underflow
    return [1 / prod(2(x̃[j] - x̃[k]) for k in 1:n if k !== j) for j in 1:n]
end

"""
	baryinterp(x, u)
	baryinterp(x, u, w)
Create a polynomial interpolant by the barycentric formula for the function values in
vector `u` at the node locations in vector `x`. If given, `w` is a vector of the
barycentric weights; otherwise it is computed from the nodes. The return value is
a callable function of the interpolation variable.
"""
function baryinterp(x, u, w=baryweights(x))
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
