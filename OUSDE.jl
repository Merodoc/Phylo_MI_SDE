using StochasticDiffEq, SciMLSensitivity, Plots
using DifferentialEquations

#Example System with OU

W = WienerProcess(0.0,0.0,0.0)
α = 0.1
σ = 0.03242599	
μ = 5.605898
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

sol