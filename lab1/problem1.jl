# CSE 5441 Numerical Methods Spring 2023
# Professor Nick Painter
#
# Author: Drew Ripberger

# Computer Problem 1.1.4
using Printf
using PrettyTables

ns = map(x -> 8.0^x, 1:10)

approx = map(x -> (1 + (1/x))^x, ns)

@printf("exp(1) = %.15f\n", exp(1))
pretty_table([ns approx], header=(["n", "approx e"]), formatters=(ft_printf("%d", [1]), ft_printf("%.15f", [2])))
