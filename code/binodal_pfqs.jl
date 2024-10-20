using Roots
using QuadGK
using NLsolve
using Plots
using LaTeXStrings

v0 = 5.0
phim = 10.0
rhom = 25.0
Dr = 1.0
lambda = 1.0

function v(rho)
    v0 * exp(-lambda * tanh((rho - rhom) / phim))
end

rf=0.25
rhos = range(0, 250, step=0.01)
plot(rhos, v.(rhos))
plot!(rhos, vstar.(rhos, rf))
plot!(rhos, v.(rhos) .* vstar.(rhos, rf) / v0)
plot!(rhos, rhos .* v.(rhos) .* vstar.(rhos, rf) / v0)


function vprime(rho)
    -lambda * v0 / phim / cosh((rho - rhom) / phim)^2 * exp(-lambda * tanh((rho - rhom) / phim))
end

function vstar(rho, rf)
    v(rho) * (1 - 0.96436447 * rho*rf^2 + 0.20893963 * rho^2*rf^4) * (1 - tanh(6.5619715 * (rho*rf^2 - 1.06965412)))/2
end

function pA(rho, rf)
    v(rho) * vstar(rho, rf) * rho / (2 * Dr)
end

function pIK(rho, rf)
    6.71615706e-2 * (exp(4.32838715 * rho*rf^2) - 1) + 5.91672633e-9 * (exp(14.8247273 * rho*rf^2) - 1)
end

function integrand(x, rf)
    -x*vprime(x)*vstar(x, rf) / (2*Dr)
end

function pG(rho, rf)
    quadgk(x -> integrand(x, rf), reference_rho, rho, rtol=1e-4)[1]
end

function p(rho, rf)
    pA(rho, rf) + pIK(rho, rf) + pG(rho, rf)
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

function Tangent(rho, rho1, rho2, F1, F2)
    F1 + (rho-rho1) * (F2-F1)/(rho2-rho1)
end


reference_rho = 25;
slope = 140
r = 0.12
rhos = range(0, 125, step=0.5);
pAs = pA.(rhos, r);
pIKs = pIK.(rhos, r);
pGs = pG.(rhos, r);
ps = pAs + pIKs + pGs;
Fraws = F.(rhos, r) ; # there must be a way to more efficiently compute this
Fs = Fraws .- rhos .* slope;

results = CommonTangentFreeEnergy(r, 5.0, 110.0)
F1 = F(results[1], r)-slope*results[1];
F2 = F(results[2], r)-slope*results[2];
println(results[3])
tangent = Tangent.(rhos, results[1], results[2], F1, F2);

plot(rhos, Fs, label=L"f(\rho)", yticks=false, dpi=300, legendfontsize=15)
plot!(rhos, tangent, label=nothing)
xlabel!(L"\rho");
# ylabel!(L"f(\rho)");
title!(L"r_f=0.12")
savefig("Free_energy_rf0.12.png")

plot(rhos[50:200], (tangent-Fs)[50:200])

plot(rhos, integrand.(rhos, r))
plot(rhos, pG.(rhos, r))
integrand(29, r)
quadgk(x -> integrand(x, r), 20, 30)

plot(rhos, ps, label="Total pressure")
plot!(rhos, pAs, label="Active pressure")
plot!(rhos, pIKs, label="Direct pressure")
plot!(rhos, pGs, label="Gradient term")
xlabel!("Density");
ylabel!("Pressure");
title!("Pressures for r_f=0.12")
savefig("Pressures_rf0.12.png")


rs = range(0.40, 0.32, step=-0.002);
number_of_rhos = length(rs);
rho1s = zeros(number_of_rhos);
rho2s = zeros(number_of_rhos);
rho1_0 = 2.0;
rho2_0 = 15.0;
for (ind, r) in enumerate(rs)
    rhos = CommonTangentFreeEnergy(r, rho1_0, rho2_0);
    rho1_0, rho2_0, norm = rhos;
    rho1s[ind] = rho1_0;
    rho2s[ind] = rho2_0;
    println("$(r) $(rho1s[ind]) $(rho2s[ind]) $(rhos[3])");
    # plot rhos
end
plot!(rho1s .* rs .* rs, 5 ./rs)
plot!(rho2s .* rs .* rs, 5 ./rs)


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
