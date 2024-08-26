using StochasticDiffEq, SciMLSensitivity, Plots
using DifferentialEquations
function lotka_volterra!(du, u, p, t)
    x, y = u
    α, β, γ, δ = p
    du[1] = dx = α * x - β * x * y
    du[2] = dy = δ * x * y - γ * y
end
function OU!(du, u, p, t)
    x = u
    α, μ, σ = p
    du[1] = α*(μ-x)
end

function stoch(du, u, p, t)
    α, μ, σ = p
    du[1]  = σ*x
end

W = WienerProcess(0.0, 0.0, 0.0)
u0 = [1.0, 1.0]
tspan = (0.0, 10.0)

function multiplicative_noise!(du, u, p, t)
    x, y = u
    du[1] = p[5] * x
    du[2] = p[6] * y
end

p = [1.5, 1.0, 3.0, 1.0, 0.3, 0.3]

prob = SDEProblem(lotka_volterra!, multiplicative_noise!, u0, tspan, p)
sol = solve(prob, SOSRI())
plot(sol)

using Statistics
ensembleprob = EnsembleProblem(prob)
@time sol = solve(ensembleprob, SOSRI(), saveat = 0.1, trajectories = 10_000)
truemean = mean(sol, dims = 3)[:,:]
truevar = var(sol, dims = 3)[:,:]

# Method of Moments 

function loss(p)
    tmp_prob = remake(prob, p=p)
    ensembleprob = EnsembleProblem(tmp_prob)
    tmp_sol = solve(ensembleprob, SOSRI(), saveat = 0.1, trajectories = 1000)
    arrsol = Array(tmp_sol)
    sum(abs2, truemean-mean(arrsol, dims = 3)) +
    0.1sum(abs2, truevar - var(arrsol, dims = 3)),
    arrsol
end

function cb2(p, l, arrsol)
    @show p, l
    means = mean(arrsol, dims = 3)[:,:]
    vars = var(arrsol, dims = 3)[:,:]
    p1 = plot(sol[1].t, means', lw = 5)
    scatter!(p1, sol[1].t, truemean')
    p2 = plot(sol[1].t, vars', lw = 5)
    scatter!(p2, sol[1].t, truevar')
    p = plot(p1, p2, layour = (2,1))
    display(p)
    false
end

using Optimization, Zygote, OptimizationOptimisers
pinit = [1.2, 0.8, 2.5, 0.8, 0.1, 0.1]
adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf,pinit)

using StochasticDiffEq, SciMLSensitivity, Plots
using DifferentialEquations

#Example System with OU

W = WienerProcess(0.0,0.0,0.0)
α = 1.0
σ = 0.1
μ = 5.0
u0 = 5.0
f(u, p, t) =  α*(μ-u)
g(u,p, t) = σ*u
prob = SDEProblem(f, g, u0, (0.0,1.0), noise = W)
sol = solve(prob, SRIW1())
plot(sol)

ensembleprob = EnsembleProblem(prob)
sol = solve(ensembleprob, EnsembleThreads(), trajectories = 1000)

using DifferentialEquations.EnsembleAnalysis
summ = EnsembleSummary(sol, 0:0.1:1)
plot(summ, labels = "Middle 95%")
summ = EnsembleSummary(sol, 0:0.1:1, quantiles = [0.25, 0.75])
plot!(summ, labels = "Middle 50%", legend = true)