# Numerical Methods Spring 2023
# Professor Nick Painter
#
# Author: Drew Ripberger

using Printf


# DividedDifference algorithm to find a_0, ... , a_n coefficients of the newton polynomial
function DividedDifference(X, Y, n)
    a = zeros(n, n)
    a[:,1] = Y[1:n]
    for i = 2:n
        for j = 1:(n-i+1)
            a[j, i] = (a[j+1,i-1]-a[j,i-1])/(X[i+j-1] - X[j])
        end
    end
    a[1,:]
end

# Helper function to evaluate a newton polynomial of the form
# p(pt) = a[0] + a[1]*(pt-X[0]) + a[2]*(pt-X[0])*(pt-X[1]) ...
function evaluate_newton(X, pt, a)
    total = 0
    n = length(a)
    temp_prod = 1
    for i = 1:n
        total += temp_prod * a[i]
        temp_prod *= pt-X[i]
    end
    total
end

# Used for Problem 4.1.15
X = [-2, -1, 0, 1, 2, 3]
Y = [1, 4, 11, 16, 13, -4]
println(DividedDifference(X, Y, 6))

# Computer Problem 4.1.3

# Since our interval is conveniently [0, 2] lets use Chebychev points
# by shifting the original domain of [-1, 1] right one

# points for the degree 10 polynomial
n = 11
X = collect(map(x -> cos(pi*((2*x + 1)/(2*n + 1))) + 1, 1:n))
Y = collect(map(x -> exp(x), X))

# interpolate polynomial using Chebychev points
p_interp = DividedDifference(X, Y, 11)


# points for comparison
n = 100
X = collect(map(x -> cos(pi*((2*x + 1)/(2*n + 1))) + 1, 1:n))
Y = collect(map(x -> exp(x), X))

p_points = collect(map(pt -> evaluate_newton(X, pt, p_interp), X))

@printf("e^x: exp(x) - p(x)\n")
for p in zip(X, Y, p_points)
    @printf("e^%.6f: %.6f - %.6f\n", p[1], p[2], p[3])
end

# Computer Problem 4.1.9

X = [1.4, 1.25]
Y = [3.7, 3.9]

# Interpolate data as a function of y

p_interp = DividedDifference(Y, X, 2)
root = evaluate_newton(Y, 0, p_interp)
@printf("root found via interpolation: x=%f\n", root)