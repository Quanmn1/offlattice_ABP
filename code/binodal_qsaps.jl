using Roots
using QuadGK
using NLsolve


# Constants and Initial Values
v0 = 5
phim = 10                               # phi scaled by rho_m
rhom = 25
lp = 5                              # Affected by Dr
lsq = 1                               # Mean square distance of the interaction kernel
# initial_guess = [2.0, 7.0]             # Initial guesses for rguess1 and rguess2
# lambda_start = 0.8                    # First value of persistence length
# lambda_stop = 1.3                     # Last value of persistence length
# lambda_step = 0.1                    # Number of l_p's to be calculated

# Define the quorum sensing function
v(rho, lambda) = v0 * exp(-lambda * tanh((rho - rhom) / phim))

# Define the function g
g(rho, lambda) = log(v(rho, lambda) * rho)

# Define κ
κ(rho, lambda) = lsq * lambda / phim / cosh((rho - rhom) / phim)^2


# Define R. R(0)=0, and R(rho>0) > 0
function R(x, lambda)
    integrand(rho) = 1 / κ(rho, lambda)
    result, _ = quadgk(integrand, rhom, x)
    return result
end

# Define rhoofR
function rhoofR(x, lambda)
    find_zero(rho -> R(rho, lambda) - x, (0,4*rhom))
end

# Define Phi
function Phi(x, lambda)
    integrand(z) = g(rhoofR(z, lambda), lambda)
    result, _ = quadgk(integrand, rhom, x)
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
    return r1, r2, result.residual_norm
end

function rho12(lambda, rguess1, rguess2)
    r1, r2 = CommonTangent(lambda, rguess1, rguess2)
    rho1 = rhoofR(r1, lambda)
    rho2 = rhoofR(r2, lambda)
    return [rho1, rho2]
end

function F(rho, lambda)
    quadgk(x -> g(x, lambda), rhom, rho)[1]
end

function CommonTangentFreeEnergy(lambda, rhoguess1, rhoguess2)
    equations(r1, r2) = [F(r1, lambda)-F(r2, lambda)-(r1-r2)*g(r1, lambda), 
                         g(r1, lambda) - g(r2, lambda)]

    # Use a solver to find roots; choose a method like NLsolve.jl or Roots.jl
    result = nlsolve(x -> equations(x[1], x[2]), [rhoguess1, rhoguess2])
    r1, r2 = result.zero
    return r1, r2, result.residual_norm
end

function Tangent(rho, rho1, rho2, F1, F2)
    F1 + (rho-rho1) * (F2-F1)/(rho2-rho1)
end

# rhos = [rho12(lambda, initial_guess...) for lambda in range(lambda_start, stop=lambda_stop, step=lambda_step)]
rho12(1.0, -300.0, 400.0)

R(7, 1.4)

rs = range(1.4, 1.5, step = 0.05)
number_of_rhos = length(rs)
rho1s = zeros(number_of_rhos)
rho2s = zeros(number_of_rhos)
rho1_0 = -4000.0
rho2_0 = 30000.0
for (ind, r) in enumerate(rs)
    rhos = CommonTangent(r, rho1_0, rho2_0);
    rho1_0, rho2_0, norm = rhos;
    rho1s[ind] = rhoofR(rho1_0, r);
    rho2s[ind] = rhoofR(rho2_0, r);
    println("$(r) $(rho1s[ind]) $(rho2s[ind]) $(rhos[3])");
    # plot rhos
end
plot(rho1s, rs)
plot!(rho2s, rs)