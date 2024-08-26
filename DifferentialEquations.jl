using DifferentialEquations

# Define the equation
f(u,p,t) = 1.01*u
# Initial condition and Timespan
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

using Plots; gr()
plot(sol)

plot(sol, linewidth=5,title="Solution to the Linear ODE with a thick line",
xaxis = "Time (t)", yaxis= "u(t) (in μM)", label = "My Thick Line!") # legend=false
plot!(sol.t, t -> 0.5*exp(1.01t), lw = 3, ls = :dash, label = "True Solution!")

# can use saveat to control the time points we check

sol = solve(prob, saveat=0.1)

#can not save start and end

sol = solve(prob, saveat = [0.2,0.7,0.9], save_start = false, save_end = false)

# or we can make it not dense 
sol = solve(prob, dense = false)

#or only the end point
sol = solve(prob, save_everystep=false)

# Lorenz Equation

function lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end

u0 = [1.0,0.0,0.0]
p = (10,28,8/3)
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!, u0, tspan,p)
sol = solve(prob)

plot(sol)

#plot x vs y vs z (Very useful)
plot(sol, vars=(1,2,3))

# DSL for Parametrized Functions

function lotka_volterra!(du,u,p,t)
    du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
    du[2] = -p[3]*u[2] + p[4]*u[1]*u[2]
end

u0 = [1.0,1.0]
p = (1.5, 1.0, 3.0, 1.0)
tspan = (0.0,10.0)
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
sol = solve(prob)
plot(sol)

# ODE by Matrix

A = [1.0 0 0 -5
    4 -2 4 -3
    -4 0 0 1
    5 -2 2 3]
u0 = rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)

sol = solve(prob)
plot(sol)

#= 
Solve the equation
du = f(u,p,t)dt + g(u,p,t)dW

where f(u,p,t) = αu
g(u,p,t) = βu 
=#

α = 1
β = 1
u₀= 1/2

f(u,p,t) = α*u
g(u,p,t) = β*u
dt = 1//2^(4)

tspan = (0.0, 1.0)
prob = SDEProblem(f,g,u₀, (0.0,1.0))

sol = solve(prob, EM(), dt = dt)
plot(sol)

# can test with an analytic Solution

f_analytic(u₀, p, t, W) = u₀ * exp((α - (β^2) / 2) * t + β * W)
ff = SDEFunction(f, g, analytic = f_analytic)
prob = SDEProblem(ff, u₀, (0.0, 1.0))

#Can solve as an ensemble problem to solve for multiple trajectories at once
ensembleprob = EnsembleProblem(prob)

sol = solve(ensembleprob,EnsembleThreads(), trajectories = 1000)

using DifferentialEquations.EnsembleAnalysis
summ = EnsembleSummary(sol, 0:0.01:1)
plot(summ, labels = "Middle 95%")
summ = EnsembleSummary(sol, 0:0.01:1, quantiles = [0.25, 0.75])
plot!(summ, labels = "Middle 50%", legend = true)

#can obtain correlations

timepoint_meancor(sol, 0.2, 0.7)

#Systems of SDEs

function lorenz(du, u,p,t)
    du[1] = 10.0(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3])-u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end

function σ_lorenz(du, u, p, t)
    du[1] = 3.0
    du[2] = 3.0
    du[3] = 3.0
end

prob_sde_lorenz = SDEProblem(lorenz, σ_lorenz, [1.0,0.0,0.0], (0.0,10.0))
sol = solve(prob_sde_lorenz)

plot(sol, idxs = (1,2,3))

# Scalar Noise

f(du, u, p, t) = (du .= u)
g(udu,u,p,t) = (du .= u)
u0 = rand(4,2)
using Plots
W = WienerProcess(0.0,0.0,0.0)
prob = SDEProblem(f,g,u0, (0.0,1.0), Noise = W)
sol = solve(prob, SRIW1())
plot(sol)