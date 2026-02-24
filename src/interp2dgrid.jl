function interp2dgrid(V, method1, method2, x, y)
    W = mapslices(v -> method1(v).(x), V; dims=[1])
    return mapslices(w -> method2(w).(y), W; dims=[2])
end
