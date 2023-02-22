# Numerical Methods Spring 2023
# Professor Nick Painter
#
# Author: Drew Ripberger

using Printf
using PrettyTables
using FFTW
using Polynomials

# Implementation of the Fast Fourier Transform and its inverse
# Note: this implementation is not well optimized and is directly
#       written from our in class derivation and psuedocode
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
    la = length(a)
    lb = length(b)
    n = max(la, lb)
    if la < n 
        a = hcat(a, zeros(n-la))
    elseif length(b) < n
        b = hcat(b, zeros(n-lb))
    end
    a = hcat(a, zeros(n))
    b = hcat(b, zeros(n))
    N = 2*n
    A = fft(a)
    B = fft(b)
    C = complex(zeros(Float64, N))
    for k = 1:N
        C[k] = A[k] * B[k]
    end
    ifft(C)
end

# Horner's method
function poly_eval(a, x)
    n = length(a)
    res = a[n]
    for k = n-1:-1:1
        res *= x
        res += a[k]
    end
    res
end

p = [10, 4, 9, 5]
q = [9, 9, 2, 5]

println(fft(p))

println(fast_polynomial_multiplication(p, q))

println(Polynomial(p)*Polynomial(q))