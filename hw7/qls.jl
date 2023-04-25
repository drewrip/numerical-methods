using LinearAlgebra

# The data we are fitting
x = [-2, -1, 0, 1, 2]
y = Dict(zip(x, [2, 1, 1, 1, 2])) # We are using this "hacky" behavior so we can treat
                                  # y like a function later when we take the inner product
yf = x -> y[x]

# Our basis functions, the monomials
g = [
     x -> 1,
     x -> x,
     x -> x^2, # You can add an arbitary number of functions to the list!
    ]

# ip is our inner product function.
# It takes two functions as input f and g,
# and it returns a function that maps a vector
# of reals to the innner product over f and g
ip(f, g) = x -> sum(map(j -> f(j)*g(j), x))

# We take the cartestian product of g with itself, ie g x g,
# in order to construct a matrix with the proper inner products
# in the right places. This is done with Iterators.product(g, g).
#
# Then we map our inner product function ip over each of the elements
# in the matrix and call the function that it returns, passing it our data x.
A = map(h -> ip(h[1], h[2])(x), Iterators.product(g, g))

# We perform a similar process with b
b = map(h -> ip(h, yf)(x), g)

# Julia allows us to use the \ operator to solve the system
c = A\b

display(c)
