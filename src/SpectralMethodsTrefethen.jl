module SpectralMethodsTrefethen

using FFTW

# Analogs for named functions in the book.
export cheb, chebfft, clencurt, gauss
include("bookfuns.jl")

# Interpolation functions for Chebyshev and Fourier discretizations.
export baryinterp, chebinterp, fourinterp, interp2d
include("interpfuns.jl")

end
