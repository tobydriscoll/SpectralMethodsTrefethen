function interp2dgrid(V, method1, method2, ξ, η)
    W = mapslices(v -> method1(v).(ξ), V; dims=[1])
    return mapslices(w -> method2(w).(η), W; dims=[2])
end
