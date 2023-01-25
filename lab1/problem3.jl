# Numerical Methods Spring 2023
# Professor Nick Painter
#
# Author: Drew Ripberger

# Computer Problem 1.2.1
using Printf
using Plots

arctan(x, n) = sum(map((w, k) -> ((-1)^(k+1))*((w^((2*k)-1))/((2*k)-1)), fill(x, n), 1:n))

x = range(-10, 10, length=50)

y1 = arctan.(x, 1)
y2 = arctan.(x, 2)
y3 = arctan.(x, 3)
y4 = arctan.(x, 4)
y5 = arctan.(x, 5)

plot(x, [y1, y2, y3, y4, y5], 
    layout=(5, 1),
    title=["1st partial sum" "2nd partial sum" "3rd partial sum" "4th partial sum" "5th partial sum"],
    size=(1000,1000),
)

xlabel!("x")
ylabel!("y")
png("partial_sums.png")