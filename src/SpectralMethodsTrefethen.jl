module SpectralMethodsTrefethen

using LinearAlgebra, FFTW

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

end
