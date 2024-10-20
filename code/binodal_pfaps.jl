using Roots
using QuadGK
using NLsolve
using Plots
using ForwardDiff
using Romberg

v0 = 5.0

# Fit with error bars
function vstar(rho)
    v0 * (1 - 0.96436447 * rho + 0.20893963 * rho^2) * (1 - tanh(6.5619715 * (rho - 1.06965412)))/2
end

function pA(rho, Dr)
    v0 * vstar(rho) * rho / (2 * Dr)
end

function pIK(rho)
    6.71615706e-2 * (exp(4.32838715 * rho) - 1) + 5.91672633e-9 * (exp(14.8247273 * rho) - 1)
end

function pIKprime(rho)
    ForwardDiff.derivative(r -> pIK(r), rho)
end

# Fit without error bars
function vstar(rho)
    v0 * exp(-rho) * exp(-2*rho^3) + rho*(1-tanh((rho-0.204)^6))^2
end

function pA(rho, Dr)
    v0 * vstar(rho) * rho / (2 * Dr)
end

function pIK(rho)
    3.65708012e-01 * (exp(2.22794865 * rho) - 1) + 2.64313115e-03 * (exp(6.59474444 * rho) - 1)
end

function p(rho, Dr)
    pA(rho, Dr) + pIK(rho)
end


function F(rho, rf)
    quadgk(x -> p(x, rf), reference_rho, rho)[1]
end

function CommonTangentFreeEnergy(rf, rhoguess1, rhoguess2)
    equations(r1, r2) = [F(r1, rf)-F(r2, rf)-(r1-r2)*p(r1, rf), 
                         p(r1, rf) - p(r2, rf)]

    # Use a solver to find roots; choose a method like NLsolve.jl or Roots.jl
    result = nlsolve(x -> equations(x[1], x[2]), [rhoguess1, rhoguess2])
    r1, r2 = result.zero
    return r1, r2, result.residual_norm
end

function Tangent(rho, rho1, rho2, F1, F2)
    F1 + (rho-rho1) * (F2-F1)/(rho2-rho1)
end
# calculate local theory and spinodal



function R(rho, Dr)
    pIK(rho)
end

function rhoofR(x, Dr)
    find_zero(y -> R(y, Dr) - x, (0.0,1.9))
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
    return r1, r2, result.residual_norm
end

function rho12(Dr, rguess1, rguess2)
    r1, r2 = CommonTangent(Dr, rguess1, rguess2)
    rho1 = rhoofR(r1, Dr)
    rho2 = rhoofR(r2, Dr)
    return [rho1, rho2]
end

# Sigma_m?


function dpdrho(rho, Dr)
    ForwardDiff.derivative(r -> p(r, Dr), rho)
end

# Find extrema by finding where the derivative is zero
function find_extrema(Dr, initial_guess1, initial_guess2)
    # Use find_zeros to find where dpdrho == 0
    critical_points = find_zeros(rho -> dpdrho(rho, Dr), initial_guess1, initial_guess2)
    return critical_points
end


# plot pressures
# D = 0.02
# Rs = range(0, 1, step=0.001);
# plot(Rs, p.(rhoofR.(Rs, D), D));
# plot(Rs, Phi.(Rs, D));

# find spinodal
# Ds = [0.05, 0.08, 0.10, 0.15]
# find_critical_points.(Ds, 0.4, 1.1)

# find binodals
reference_rho = 0
Drs = range(0.02, 0.30, step = 0.005)
number_of_rhos = length(Drs)
rho1s = zeros(number_of_rhos);
rho2s = zeros(number_of_rhos);
ps = zeros(number_of_rhos)
r1_0 = 0.01
r2_0 = 60.0
for (ind, Dr) in enumerate(Drs)
    rs = CommonTangent(Dr, r1_0, r2_0);
    r1_0, r2_0, norm = rs;
    rho1s[ind] = rhoofR(r1_0, Dr);
    rho2s[ind] = rhoofR(r2_0, Dr);
    ps[ind] = p(rho1s[ind], Dr);
    println("$(Dr) $(rho1s[ind]) $(rho2s[ind]) $(rs[3])");
    # plot rhos
end
plot(rho1s, 5 ./Drs)
plot!(rho2s, 5 ./Drs)

include("RungeKutta.jl")

# Influence of Sigma_m:
# Use binodals. Solve for d\rho/dx.
# Integrate to calculate correction to P_{co}, hence correction to binodals
# Iterate if necessary

# Take a particular Dr, with that a particular P_{co} to zeroth order.
function a(rho)
    vstar(rho) * pIKprime(rho)
end

function b(rho, Dr)
    16/3 * (p(rho, Dr) - p_co)/(vstar(rho) * (v0/Dr)^2)
end

# construct the system of ODEs
function forward(t, vars, params)
    Dr = params[1]
    rho = vars[1]
    lambda = vars[2]
    return [lambda / a(rho), b(rho, Dr)]
end

# stop condition: the change in rho becomes small
function stop_condition(t, vars)
    rho = vars[1]
    lambda = vars[2]
    return abs(lambda) < 0.001
end

# rho_profile[i]: 1 is rho, 2 is a(rho) * drho/dx
function profile_integrand(rho_profile)
    rho = rho_profile[1]
    rho_grad = rho_profile[2] / a(rho)
    return rho_grad^3 * pIKprime(rho)^2 * (pIKprime(rho)/rho - pIK(rho)/rho^2)
end

function maxwell(p_co, Dr)
    rho_gas_new = find_zero(x -> p(x, Dr) - p_co, rho_gas)
    rho_liquid_new = find_zero(x -> p(x, Dr) - p_co, rho_liquid)
    return quadgk(x -> (p(x, Dr) - p_co)*pIKprime(x), rho_gas_new, rho_liquid_new)[1]
end

rhos = range(0, 1.4, 1000)
rho1s_new = zeros(number_of_rhos);
rho2s_new = zeros(number_of_rhos);

# start the ODE close to the bulk but not in the bulk
for (ind, Dr) in enumerate(Drs)
    # ind = 1
    # Dr = Drs[ind]
    # relevant parameters
    p_co = ps[ind]
    rho_gas = rho1s[ind]
    rho_liquid = rho2s[ind]
    # do the ODE. remember to try different initial conditions!
    xs, rho_profiles = rk4(forward, [rho_gas+0.001,0.001], 0.0, stop_condition, params=[Dr], step_size=0.001)
    xs = LinRange(xs[1], xs[end], length(xs))
    integral, error = romberg(xs, profile_integrand.(rho_profiles))

    # see how the result changes when turning xi up
    xi = 100
    int = integral / 2 / v0 / Dr * xi
    # println(int)
    # find the modified gas and liquid density
    p_co_new = find_zero(x -> maxwell(x, Dr) - int, p_co)
    rho1s_new[ind] = find_zero(x -> p(x, Dr) - p_co_new, rho_gas)
    rho2s_new[ind] = find_zero(x -> p(x, Dr) - p_co_new, rho_liquid)
    println("$(Dr) $(rho1s_new[ind]) $(rho2s_new[ind])");
end

plot!(rho1s_new, 5 ./Drs)
plot!(rho2s_new, 5 ./Drs)

