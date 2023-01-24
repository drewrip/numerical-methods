# CSE 5441 Numerical Methods Spring 2023
# Professor Nick Painter
#
# Author: Drew Ripberger

# Computer Problem 1.2.1
using Printf
using Plots

arctan(x, n) = sum(
    map((w, k) -> ((-1)^(k+1))*((w^((2*k)-1))/((2*k)-1)), fill(x, n), 1:n))

x = range(-10, 10, length=25)
y = arctan.(x, 5)

plot(x, y)
xlabel!("x")
ylabel!("y")
png("prob3.png")