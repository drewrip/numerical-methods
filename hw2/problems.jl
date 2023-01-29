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