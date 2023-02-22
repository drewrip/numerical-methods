# Numerical Methods Spring 2023
# Professor Nick Painter
#
# Author: Drew Ripberger

using Printf
using PrettyTables
using FFTW
using Polynomials

function fft(a)
    n = length(a)
    if n == 1
        a
    else
        a_e = a[1:2:end]
        a_o = a[2:2:end]
        y_o = fft(a_o)
        y_e = fft(a_e)

        w = 1
        psi = exp((2 * pi * 1im)/n)
        y = complex(zeros(Float64, n))
        for k = 1:(floor(Int64, n/2))
            println(k)
            y[k] = y_e[k] + w*y_o[k]
            y[k + floor(Int64, n/2)] = y_e[k] - w*y_o[k]
            w = w * psi
        end
        y
    end
end

function ifft_helper(a)
    n = length(a)
    if n == 1
        a
    else
        a_e = a[1:2:end]
        a_o = a[2:2:end]

        y_o = ifft_helper(a_o)
        y_e = ifft_helper(a_e)

        w = 1
        psi = exp((-2 * pi * 1im)/n)
        y = complex(zeros(Float64, n))
        for k = 1:(floor(Int64, n/2))
            y[k] = y_e[k] + w*y_o[k]
            y[k + floor(Int64, n/2)] = y_e[k] - w*y_o[k]
            w = w * psi
        end
        y
    end
end

function ifft(a)
    n = length(a)
    y = ifft_helper(a)
    y/n
end

function fast_polynomial_multiplication(a, b)
    # Pad the input polynomials if needed in future
    n = length(a)
    A = fft(a)
    B = fft(b)
    C = complex(zeros(Float64, n))
    for k = 1:n
        C[k] = A[k] * B[k]
    end
    ifft(C)
end

function poly_eval(a, x)
    res = 0+0im
    n = length(a)
    xp = 1
    for k = 1:n
        res += xp*a[k]
        xp *= x
    end
    res
end

p = [10, 4, 9, 5]
q = [9, 9, 2, 5]

println(fast_polynomial_multiplication(p, q))

println(Polynomial(p)*Polynomial(q))