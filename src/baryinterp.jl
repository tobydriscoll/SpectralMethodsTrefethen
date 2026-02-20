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
