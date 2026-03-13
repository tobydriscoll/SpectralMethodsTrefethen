using Test, LinearAlgebra, SpectralMethodsTrefethen

@testset "fourier" begin
    for N in [26, 32]
        x, D, D2 = fourier(N)
        @test x[end] ≈ 2π
        v = @. exp(cos(x))
        @test D * v ≈ -sin.(x) .* v
        @test D2 * v ≈ (sin.(x).^2 - cos.(x)) .* v
        v = @. exp(cis(x))
        @test D * v ≈ im * cis.(x) .* v
    end
end

@testset "fourinterp" begin
    ξ = range(0, 2π, 23)
    for N in [4, 8]
        x, _ = fourier(N)
        u = fourinterp(sin.(x))
        @test u.(ξ) ≈ sin.(ξ) atol=1e-8
    end
    for N in [28, 44]
        x, _ = fourier(N)
        for f in (x -> exp(2cos(x)), x -> cis.(-x))
            u = fourinterp(f.(x))
            @test u.(ξ) ≈ f.(ξ) atol=1e-8
        end
    end
end

@testset "cheb" begin
    for N in 0:3
        D, x = cheb(N)
        @test size(D) == (N+1, N+1)
        @test length(x) == N+1
        if N > 0
            @test D * x ≈ ones(N+1)
        else
            @test D * x ≈ [0]
        end
    end
    for N in [17, 22, 39]
        D, x = cheb(N)
        @test x[1] ≈ 1
        @test D * cos.(x) ≈ -sin.(x)
        @test D * cis.(2x) ≈ 2im * cis.(2x)
    end
end

@testset "chebfft" begin
    @test chebfft([1.]) ≈ [0]
    for N in 1:3
        D, x = cheb(N)
        @test chebfft(x) ≈ ones(N+1)
    end
    for N in [17, 22, 31]
        D, x = cheb(N)
        @test chebfft(cos.(x)) ≈ -sin.(x)
        @test chebfft(cis.(2x)) ≈ 2im * cis.(2x)
    end
end

@testset "chebinterp" begin
    ξ = -1:0.1:1
    for N in 1:3
        D, x = cheb(N)
        u = chebinterp(x)
        @test u.(ξ) ≈ ξ
    end
    for N in [17, 22, 31]
        D, x = cheb(N)
        for f in (cos, x -> cis.(2x))
            v = f.(x)
            u = chebinterp(v)
            @test u.(ξ) ≈ f.(ξ)
        end
    end
end

@testset "integration by $method" for method in (clencurt, gauss)
    for N in 2:5
        x, w = method(N)
        @test sum(w) ≈ 2
        @test sum(w .* x) ≈ 0 atol=1e-8
        @test sum(w .* x.^2) ≈ 2/3
    end
    for N in [15, 26]
        x, w = method(N)
        for (f, I) in [(x -> exp(-4x), sinh(4) / 2), (sech, 2atan(sinh(1))) ]
            @test sum(w .* f.(x)) ≈ I atol=1e-8
        end
    end
end

@testset "2-D interp" begin
    f(x, y) = exp(-x) + tanh(x + y)
    _, x = cheb(16)
    _, y = cheb(21)
    V = [f(x, y) for x in x, y in y]
    ξ = η = -1:0.025:1
    @test interp2dgrid(V, chebinterp, chebinterp, ξ, η) ≈ [f(x, y) for x in ξ, y in η]
end