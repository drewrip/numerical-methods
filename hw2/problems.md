# HW 2 (CSE 5361 Painter) - Drew Ripberger
The solutions are written in Julia.

## Code

```julia
# Numerical Methods Spring 2023
# Professor Nick Painter
#
# Author: Drew Ripberger

using Printf
using PrettyTables
using Polynomials

# Computer Problem 3.1.1

# Find the roots of f
function bisection(f, lower, upper, N, epsilon)
    guess = (upper-lower)/2 + lower
    f_guess = f(guess)
    if abs(f_guess) < epsilon || N < 1
        guess
    elseif sign(f_guess) == sign(f(lower))
        bisection(f, guess, upper, N-1, epsilon)
    elseif sign(f_guess) == sign(f(upper))
        bisection(f, lower, guess, N-1, epsilon)
    end
end

f = (x) -> x^3 - x^2 - 2*x + 1

a = -100000
b = 100000

println("root of f(x) at x=", bisection(f, a, b, 100, .01))

# Computer Problem 3.1.4

g = (x) -> 6 + 3*(x^2) + 2*(x^3) - 6*(exp(x)-x)

x = -1
y = 1

println("root of g(x) at x=", bisection(g, x, y, 100, .01))

# Computer Problem 3.1.10

h = Polynomial([40320, -109584, 118124, -67284, 22449, -4536, 546, -36, 1])
println(h)
println("roots of the polynomial:")
println(roots(h))
h_prime = Polynomial([40320, -109584, 118124, -67284, 22449, -4536, 546, -37, 1])
println(h_prime)
println("roots of the new polynomial:")
println(roots(h_prime))
```

## Output

\small
```
root of f(x) at x=-1.245737075805664
root of g(x) at x=0.0
40320 - 109584*x + 118124*x^2 - 67284*x^3 + 22449*x^4 - 4536*x^5 + 546*x^6 - 36*x^7 + x^8
roots of the polynomial:
[0.9999999999999909, 1.9999999999999936, 3.000000000000528, 3.9999999999987774, 4.999999999999056,
    6.000000000006543, 6.999999999991528, 8.00000000000362]
40320 - 109584*x + 118124*x^2 - 67284*x^3 + 22449*x^4 - 4536*x^5 + 546*x^6 - 37*x^7 + x^8
roots of the new polynomial:
ComplexF64[0.999801963896044 + 0.0im, 2.0843875381075616 - 0.24935240473949752im, 
    2.0843875381075616 + 0.24935240473949752im, 2.821038133239228 - 1.7281215850060632im,
    2.821038133239228 + 1.7281215850060632im, 5.035095810228782 - 5.149749378225431im,
    5.035095810228782 + 5.149749378225431im, 16.119155072952807 + 0.0im]
```
\normalsize