using Roots
using QuadGK
using NLsolve
using Plots
using ForwardDiff

v0 = 5.0

function vstar(rho)
    v0 * (1 - 0.95902466 * rho + 0.21347265 * rho^2) * (1 - tanh(5.96025284 * (rho - 1.04069886)))/2
end

function pA(rho, Dr)
    v0 * vstar(rho) * rho / (2 * Dr)
end

function pIK(rho)
    6.42305094e-2 * (exp(4.35125642 * rho) - 1) + 6.94245881e-9 * (exp(14.7539197 * rho) - 1)
end

function p(rho, Dr)
    pA(rho, Dr) + pIK(rho)
end

function R(rho, Dr)
    pIK(rho)
end

function rhoofR(x, Dr)
    find_zero(y -> R(y, Dr) - x, (0.0,2.9))
end

function Phi(x, Dr)
    quadgk(z -> p(rhoofR(z, Dr), Dr), 1, x)[1]
end

function q(r, r1, Dr)
    Phi(r1, Dr) + (r - r1) * p(rhoofR(r1, Dr), Dr)
end

function CommonTangent(Dr, rguess1, rguess2)
    # Define the system of equations
    equations(r1, r2) = [q(r2, r1, Dr) - Phi(r2, Dr), 
                         p(rhoofR(r1, Dr), Dr) - p(rhoofR(r2, Dr), Dr)]

    # Use a solver to find roots; choose a method like NLsolve.jl or Roots.jl
    result = nlsolve(x -> equations(x[1], x[2]), [rguess1, rguess2])

    r1, r2 = result.zero
    return r1, r2
end

function rho12(Dr, rguess1, rguess2)
    r1, r2 = CommonTangent(Dr, rguess1, rguess2)
    rho1 = rhoofR(r1, Dr)
    rho2 = rhoofR(r2, Dr)
    return [rho1, rho2]
end

function dpdrho(rho, Dr)
    ForwardDiff.derivative(r -> p(r, Dr), rho)
end

# Find critical points by finding where the derivative is zero
function find_critical_points(Dr, initial_guess1, initial_guess2)
    # Use find_zeros to find where dpdrho == 0
    critical_points = find_zeros(rho -> dpdrho(rho, Dr), initial_guess1, initial_guess2)
    return critical_points
end



D = 0.01
rhos = range(0,1.4,step=0.01)
Rs = range(0, 75, step=0.1)
plot(Rs, p.(rhoofR.(Rs, D), D))
plot(Rs, Phi.(Rs, D))

Ds = [0.05, 0.08, 0.10, 0.15]
find_critical_points.(Ds, 0.4, 1.1)

plot(rhos, p.(rhos, 0.05))

Ds = [0.20, 0.25, 0.30]
Ds=[0.01,0.02]
Ds=[0.01]
results = rho12.(Ds, 0.1, 200.0)
for (i, res) in enumerate(results)
    println("D = $(Ds[i]): rho12 = $res")
end
