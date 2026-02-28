module SpectralMethodsTrefethen

using LinearAlgebra, FFTW, CairoMakie, ForwardDiff, LaTeXStrings, Polynomials, Printf, SpecialFunctions, ToeplitzMatrices

# Analogs for named functions in the book.
export cheb, chebfft, clencurt, gauss
include("cheb.jl")
include("chebfft.jl")
include("clencurt.jl")
include("gauss.jl")

# Interpolation functions for Chebyshev and Fourier discretizations.
export baryinterp, chebinterp, fourinterp, interp2dgrid
include("baryinterp.jl")
include("chebinterp.jl")
include("fourinterp.jl")
include("interp2dgrid.jl")

# Main programs from the book.
export p1, p2, p3, p4, p5, p6, p6anim, p7, p8, p9, p10
export p11, p12, p13, p14, p15, p16, p17, p18, p19, p19anim, p20anim
export p21, p22, p23, p24, p25, p26, p27, p27anim, p28, p29, p30
export p31, p32, p33, p34, p34anim, p35, p35anim, p36, p37anim, p38, p39, p40, p40anim
for file in filter(endswith("jl"), readdir(joinpath(@__DIR__, "scripts")))
    include(joinpath("scripts", file))
    precompile(eval(Symbol(splitext(file)[1])), ())
end

end    # module