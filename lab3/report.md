# Lab 3 (CSE 5361 Painter) - Drew Ripberger

The solutions are written in Julia.

## Code

```julia
using Printf
using LinearAlgebra


function factorize_algo(A, L_diag, U_diag)
    n = size(A, 1)
    L = zeros((n,n))
    U = zeros((n,n))
    for i in 1:n
        diag_entry = (A[i, i] - sum(
            map(x -> x[1]*x[2], zip(L[i, begin:i-1], U[begin:i-1, i])),
                init=0
            )
        )
        if L_diag[i] == 0
            L[i, i] = diag_entry/U_diag[i]
            U[i, i] = U_diag[i]
        else
            L[i, i] = L_diag[i]
            U[i, i] = diag_entry/L_diag[i]
        end

        for j in i:n
            L[j, i] = (1/U[i, i])*(A[j, i] - sum(
                map(x -> x[1]*x[2], zip(L[j, begin:i-1], U[begin:i-1, i])),
                    init=0
                )
            )
            U[i, j] = (1/L[i, i])*(A[i, j] - sum(
                map(x -> x[1]*x[2], zip(L[i, begin:i-1], U[begin:i-1, j])),
                    init=0
                )
            )
        end
    end
    L, U   
end

# Diagonal of L is fixed to 1's
function doolittle(A)
    factorize_algo(A, ones(size(A, 1)), zeros(size(A, 1)))
end

# Diagonal of U is fixed to 1's
function crout(A)
    factorize_algo(A, zeros(size(A, 1)), ones(size(A, 1)))    
end

# This is crazy, this is basically identical to our psuedocode LOL
function eig_power(A, x, psi, kmax)
    r = 0
    y = zeros(size(A, 1))
    for k in 1:kmax
        y = A*x
        r = psi(y)/psi(x)
        x = y/norm(y, 2)
    end
    r, x 
end

# Solve system in matrix A, assuming it is lower triangular
function forwardsub(A, b)
    n = size(A, 1)
    known = [b[1]/A[1, 1]]
    for i in 2:n
        x = (b[i] - dot(known, A[i, 1:i-1]))/A[i, i] 
        push!(known, x)
    end
    known
end

# Solve system in matrix A, assuming it is upper triangular
function backsub(A, b)
    n = size(A, 1)
    known = [b[n]/A[n, n]]
    for i in n-1:-1:1
        x = (b[i] - dot(known, A[i, i+1:n]))/A[i, i]
        push!(known, x)
    end
    reverse(known)
end

function mv_newtons(f, Df, x0)
    x = x0
    for i in 1:100
        L, U = doolittle(map(g -> g(x), Df))
        b = map(g -> -g(x), f)
        y = forwardsub(L, b)
        hk = forwardsub(U, y)
        x = x + hk
    end
    x
end

M = [5 7 6 5; 7 10 8 7; 6 8 10 9; 5 7 9 10]

L, U = doolittle(M)

println("Doolittle factorization:")
println("L =")
display(L)
println("U =")
display(U)
println("L*U =")
display(L * U)

L, U = crout(M)

println("Crout factorization:")
println("L =")
display(L)
println("U =")
display(U)
println("L*U =")
display(L * U)

L, U = factorize_algo(M, [1 0 1 0], [0 1 0 1])

println("[1 ? 1 ?], [? 1 ? 1] factorization:")
println("L =")
display(L)
println("U =")
display(U)
println("L*U =")
display(L * U)

L, U = factorize_algo(M, [0 1 0 1], [1 0 1 0])

println("[? 1 ? 1], [1 ? 1 ?] factorization:")
println("L =")
display(L)
println("U =")
display(U)
println("L*U =")
display(L * U)

L, U = factorize_algo(M, [0 0 7 9], [3 5 0 0])

println("[? ? 7 9], [3 5 ? ?] factorization:")
println("L =")
display(L)
println("U =")
display(U)
println("L*U =")
display(L*U)

println("========================================")

r, x = eig_power(M, [1.0; 2.0; 3.0; 4.0], (x) -> x[3], 100)
println("dominant eigenvalue of M = ", r)
println("corresponding dominant eigenvector = ")
display(x)

println("mv newtons method, x = ")
x = mv_newtons(
        [
         (x) -> x[1]^2 + x[2]^2 - 5,
         (x) -> x[1] + 2*x[2] - 4,
        ],
        [
         ((x) -> 2*x[1]) ((x) -> 2*x[2]);
         ((x) -> 1) ((x) -> 2);
        ],
        [1, 1],
    )
display(x)
```

## Output

```
Doolittle factorization:
L =
4×4 Matrix{Float64}:
 1.0   0.0  0.0  0.0
 1.4   1.0  0.0  0.0
 1.2  -2.0  1.0  0.0
 1.0   0.0  1.5  1.0
U =
4×4 Matrix{Float64}:
 5.0  7.0   6.0   5.0
 0.0  0.2  -0.4  -8.88178e-16
 0.0  0.0   2.0   3.0
 0.0  0.0   0.0   0.5
L*U =
4×4 Matrix{Float64}:
 5.0   7.0   6.0   5.0
 7.0  10.0   8.0   7.0
 6.0   8.0  10.0   9.0
 5.0   7.0   9.0  10.0
Crout factorization:
L =
4×4 Matrix{Float64}:
 5.0   0.0          0.0  0.0
 7.0   0.2          0.0  0.0
 6.0  -0.4          2.0  0.0
 5.0  -8.88178e-16  3.0  0.5
U =
4×4 Matrix{Float64}:
 1.0  1.4   1.2  1.0
 0.0  1.0  -2.0  0.0
 0.0  0.0   1.0  1.5
 0.0  0.0   0.0  1.0
L*U =
4×4 Matrix{Float64}:
 5.0   7.0   6.0   5.0
 7.0  10.0   8.0   7.0
 6.0   8.0  10.0   9.0
 5.0   7.0   9.0  10.0
[1 ? 1 ?], [? 1 ? 1] factorization:
L =
4×4 Matrix{Float64}:
 1.0   0.0  0.0  0.0
 1.4   0.2  0.0  0.0
 1.2  -0.4  1.0  0.0
 1.0   0.0  1.5  0.5
U =
4×4 Matrix{Float64}:
 5.0  7.0   6.0   5.0
 0.0  1.0  -2.0  -4.44089e-15
 0.0  0.0   2.0   3.0
 0.0  0.0   0.0   1.0
L*U =
4×4 Matrix{Float64}:
 5.0   7.0   6.0   5.0
 7.0  10.0   8.0   7.0
 6.0   8.0  10.0   9.0
 5.0   7.0   9.0  10.0
[? 1 ? 1], [1 ? 1 ?] factorization:
L =
4×4 Matrix{Float64}:
 5.0   0.0          0.0  0.0
 7.0   1.0          0.0  0.0
 6.0  -2.0          2.0  0.0
 5.0  -4.44089e-15  3.0  1.0
U =
4×4 Matrix{Float64}:
 1.0  1.4   1.2  1.0
 0.0  0.2  -0.4  0.0
 0.0  0.0   1.0  1.5
 0.0  0.0   0.0  0.5
L*U =
4×4 Matrix{Float64}:
 5.0   7.0   6.0   5.0
 7.0  10.0   8.0   7.0
 6.0   8.0  10.0   9.0
 5.0   7.0   9.0  10.0
[? ? 7 9], [3 5 ? ?] factorization:
L =
4×4 Matrix{Float64}:
 1.66667   0.0    0.0  0.0
 2.33333   0.04   0.0  0.0
 2.0      -0.08   7.0  0.0
 1.66667   0.0   10.5  9.0
U =
4×4 Matrix{Float64}:
 3.0  4.2    3.6       3.0
 0.0  5.0  -10.0       0.0
 0.0  0.0    0.285714  0.428571
 0.0  0.0    0.0       0.0555556
L*U =
4×4 Matrix{Float64}:
 5.0   7.0   6.0   5.0
 7.0  10.0   8.0   7.0
 6.0   8.0  10.0   9.0
 5.0   7.0   9.0  10.0
========================================
dominant eigenvalue of M = 30.288685345802126
corresponding dominant eigenvector =
4-element Vector{Float64}:
 0.38026207439071347
 0.5285678495286417
 0.5519548496316626
 0.5209247807436572
mv newtons method, x =
2-element Vector{Float64}:
 2.0
 1.0
```
