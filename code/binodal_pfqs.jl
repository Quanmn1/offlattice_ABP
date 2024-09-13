using Roots
using QuadGK
using NLsolve
using Plots

v0 = 5.0
phim = 10.0
rhom = 25.0
Dr = 1.0
lambda = 1.0

function v(rho)
    v0 * exp(-lambda * tanh((rho - rhom) / phim))
end

function vprime(rho)
    -lambda * v0 / phim / cosh((rho - rhom) / phim)^2 * exp(-lambda * tanh((rho - rhom) / phim))
end

function vstar(rho, rf)
    v(rho) * (1 - 0.95902466 * rho*rf^2 + 0.21347265 * rho^2*rf^4) * (1 - tanh(5.96025284 * (rho*rf^2 - 1.04069886)))/2
end

function pA(rho, rf)
    v(rho) * vstar(rho, rf) * rho / (2 * Dr)
end

function pIK(rho, rf)
    6.42305094e-2 * (exp(4.35125642 * rho*rf^2) - 1) + 6.94245881e-9 * (exp(14.7539197 * rho*rf^2) - 1)
end

function pG(rho, rf)
    quadgk(x -> -x*vprime(x)*vstar(x, rf) / (2*Dr), 100, rho)[1]
end

function p(rho, rf)
    pA(rho, rf) + pIK(rho, rf) + pG(rho, rf)
end

function F(rho, rf)
    quadgk(x -> p(x, rf), 100, rho)[1]
end

function CommonTangentFreeEnergy(rf, rhoguess1, rhoguess2)
    equations(r1, r2) = [F(r1, rf)-F(r2, rf)-(r1-r2)*p(r1, rf), 
                         p(r1, rf) - p(r2, rf)]

    # Use a solver to find roots; choose a method like NLsolve.jl or Roots.jl
    result = nlsolve(x -> equations(x[1], x[2]), [rhoguess1, rhoguess2])

    r1, r2 = result.zero
    return r1, r2
end

r = 0.07
rhos = range(0, 120, step=0.1)
plot(rhos, p.(rhos, r), label="Total pressure")
plot!(rhos, pA.(rhos, r), label="Active pressure")
plot!(rhos, pIK.(rhos, r), label="Direct pressure")
plot!(rhos, pG.(rhos, r), label="Gradient term")
xlabel!("Density")
ylabel!("Pressure")
title!("Pressures for r_f=0.12")
savefig("Pressures_rf0.12.png")

plot(rhos, F.(rhos, r), label="F")
xlabel!("Density")
ylabel!("Free energy")
title!("Free energy for r_f=0.12")
savefig("Free_energy_rf0.12.png")


CommonTangentFreeEnergy(r, 10.0, 250.0)

rs = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18]
rs = [0.05]
results = CommonTangentFreeEnergy.(rs, 10.0, 90.0)
for (i, res) in enumerate(results)
    println("r_f = $(rs[i]): rho12 = $res")
end

# function R(rho, Dr)
#     pIK(rho)
# end

function rhoofR(x, Dr)
    find_zero(y -> R(y, Dr) - x, (0.01,1.5))
end

function Phi(x, Dr)
    quadgk(z -> p(rhoofR(z, Dr), Dr), 1, x)[1]
end

function q(r, r1, Dr)
    Phi(r1, Dr) + (r - r1) * p(rhoofR(r1, Dr), Dr)
end

function CommonTangent(rf, rguess1, rguess2)
    # Define the system of equations
    equations(r1, r2) = [q(r2, r1, rf) - Phi(r2, rf), 
                         p(rhoofR(r1, rf), rf) - p(rhoofR(r2, rf), rf)]

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


rho12(0.01, 0.1, 60.0)
