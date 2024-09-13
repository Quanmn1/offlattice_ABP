using Roots
using QuadGK
using NLsolve


# Constants and Initial Values
phim = 2/5                               # phi scaled by rho_m
ε = 0                                 # Epsilon (offset of the QS kernel)
lp = 5                              # Affected by Dr
lsq = 1                               # Mean square distance of the interaction kernel
initial_guess = [2.0, 7.0]             # Initial guesses for rguess1 and rguess2
lambda_start = 0.8                    # First value of persistence length
lambda_stop = 1.3                     # Last value of persistence length
lambda_step = 0.1                    # Number of l_p's to be calculated

# Define the quorum sensing function
v(rho, lambda) = exp(-lambda * tanh((rho - 1) / phim))

# Define the function g
g(rho, lambda) = log(v(rho, lambda) * rho)

# Define κ
κ(rho, lambda) = lsq * lambda / phim / cosh((rho - 1) / phim)^2



# Define R. R(0)=0, and R(rho>0) > 0
function R(x, lambda)
    integrand(rho) = 1 / κ(rho, lambda)
    result, _ = quadgk(integrand, 0, x)
    return result
end

# Define rhoofR
function rhoofR(x, lambda)
    find_zero(rho -> R(rho, lambda) - x, (0,2.3))
end

# Define Phi
function Phi(x, lambda)
    integrand(z) = g(rhoofR(z, lambda), lambda)
    result, _ = quadgk(integrand, 1, x)
    return result
end

# Define q
q(R, r1, lambda) = Phi(r1, lambda) + (R - r1) * g(rhoofR(r1, lambda), lambda)

function CommonTangent(lambda, rguess1, rguess2)
    # Define the system of equations
    equations(r1, r2) = [q(r2, r1, lambda) - Phi(r2, lambda), 
                         g(rhoofR(r1, lambda), lambda) - g(rhoofR(r2, lambda), lambda)]

    # Use a solver to find roots; choose a method like NLsolve.jl or Roots.jl
    result = nlsolve(x -> equations(x[1], x[2]), [rguess1, rguess2])

    r1, r2 = result.zero
    return r1, r2
end

function rho12(lambda, rguess1, rguess2)
    r1, r2 = CommonTangent(lambda, rguess1, rguess2)
    rho1 = rhoofR(r1, lambda)
    rho2 = rhoofR(r2, lambda)
    return [rho1, rho2]
end

rhos = [rho12(lambda, initial_guess...) for lambda in range(lambda_start, stop=lambda_stop, step=lambda_step)]
# rho12(1.3, initial_guess...)