# Numerical Methods Spring 2023
# Professor Nick Painter
#
# Author: Drew Ripberger

# Computer Problem 1.2.1
using Printf

function quad(a, b, c)
    disc = sqrt((b*b) - (4*a*c))
    ((-b-disc)/(2*a), (-b+disc)/(2*a))
end

println("x^2 + (10^8)x + 1 = 0 -> x=", quad(1, 100000000, 1))
println("x^2 + (10^8)x + 10^8 = 0 -> x=", quad(1, 100000000, 100000000))