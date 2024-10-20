using Roots
using QuadGK
using NLsolve
using Plots

v0 = 5.0
rhom = 25.0
Dr = 1.0
lambda = 1.0
reference_rho = 24.5
rf = 0

# condensation
phim = 1
function v(rho)
    @. if rho < rhom - phim
        return v0
    elseif rho < rhom
        return v0 * (rhom - rho) / phim
    else
        return 0
    end
end

function vprime(rho)
    @. if rhom - phim <= rho <= rhom
        return - v0/phim
    else
        return 0
    end

end

# start with QSAPs only
function integrand(x)
    -x*vprime(x)*v(x) / (2*Dr)
end

function pA(rho, rf)
    rho*v(rho)^2 / (2*Dr) + quadgk(x -> integrand(x), reference_rho, rho, rtol=1e-4)[1]
end

function p(rho, rf)
    pA(rho, rf)
end

function F(rho, rf)
    quadgk(x -> p(x, rf), reference_rho, rho)[1]
end

function CommonTangentFreeEnergy(rf, rhoguess1, rhoguess2)
    equations(r1, r2) = [F(r1, rf)-F(r2, rf)-(r1-r2)*p(r1, rf), 
                         p(r1, rf) - p(r2, rf)]

    # Use a solver to find roots; choose a method like NLsolve.jl or Roots.jl
    result = nlsolve(x -> equations(x[1], x[2]), [rhoguess1, rhoguess2])
    # println(result)
    r1, r2 = result.zero
    return r1, r2, result.residual_norm
end

rhos = range(0, 60, step=0.1)
plot(rhos, p.(rhos, rf))
plot(rhos, F.(rhos, rf))