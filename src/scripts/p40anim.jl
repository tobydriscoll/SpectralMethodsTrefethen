using CairoMakie, Printf, LaTeXStrings
using LinearAlgebra, SpectralMethodsTrefethen
"p40 - eigenvalues of Orr–Sommerfeld operator (compare p38)"
function p40anim(Re=5772)
    function oseigs(N)
        # 2nd- and 4th-order differentiation matrices:
        D, x = cheb(N)
        D² = (D^2)[2:N, 2:N]
        S = Diagonal([0; (@. 1 / (1 - x[2:N]^2)); 0])
        D⁴ = (Diagonal(@. 1 - x^2) * D^4 - Diagonal(8x) * D^3 - 12D^2) * S
        D⁴ = D⁴[2:N, 2:N]

        # Orr–Sommerfeld operators A,B and generalized eigenvalues:
        B = D² - I
        A = (D⁴ - 2D² + I) / Re - 2im * I - 1im * Diagonal(@. 1 - x[2:N]^2) * B
        return eigvals(A, B)
    end
    fig = Figure(size=(500, 500))
    N = Observable(10)
    λ = @lift oseigs($N)
    Reλ = @lift real($λ)
    Imλ = @lift imag($λ)
    str = @lift @sprintf("\$N = %d,\\; \\lambda_\\mathrm{max}\$ = %.11f", $N, maximum($Reλ))
    title = @lift latexstring($str)
    ax = Axis(fig[1, 1]; aspect=1, limits=(-0.8, 0.2, -1, 0), title)
    vlines!(ax, [0], color=:gray)
    scatter!(ax, Reλ, Imλ; color=Reλ, colormap=:Blues, colorrange=(-0.8, 0) )
    anim = record(fig, "p40anim.mp4", 20:4:100; framerate=2, loop=-1) do n
        N[] = n
    end
    return anim
end