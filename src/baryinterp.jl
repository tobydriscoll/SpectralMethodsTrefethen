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
	baryinterp(x, v)
	baryinterp(x, v, w)
Create a polynomial interpolant by the barycentric formula for the function values in
vector `v` at the node locations in vector `x`. If given, `w` is a vector of the
barycentric weights; otherwise, it is computed from the nodes. The return value is
a callable function of the interpolation variable.
"""
function baryinterp(x, v, w=baryweights(x))
    return function(ξ)
        numer = denom = 0.0
        for i in eachindex(x)
            term = w[i] / (ξ - x[i])
            if isinf(term) || isnan(term)    # did we just divide by zero?
                return v[i]                  # return value at node
            end
            denom += term
            numer = muladd(term, v[i], numer)
        end
        return numer / denom
    end
end