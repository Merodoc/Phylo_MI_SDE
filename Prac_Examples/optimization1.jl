using OrdinaryDiffEq, 
Optimization, OptimizationPolyalgorithms, SciMLSensitivity,
Zygote, Plots


function lotka_volterra!(du,u,p,t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α * x - β * x * y
    du[2] = dy = -δ * y + γ * x * y
end

# Initial Condition
u0 = [1.0, 1.0]

# Simulation interval and intermediary points
tspan = (0.0, 10.0)
tsteps = (0.0:0.1:10.0)

# LV equation parameter. p 
p = [1.5, 1.0, 3.0, 1.0]

prob = ODEProblem(lotka_volterra!, u0, tspan, p)
sol = solve(prob, Tsit5())

#Plot the Solution

using Plots
plot(sol)
savefig("LV_ode.png")

function loss(p)
    sol = solve(prob, Tsit5(), p = p, saveat = tsteps)
    loss = sum(abs2, sol .- 1)
    return loss, sol
end

callback = function(p, l, pred)
    display(l)
    plt = plot(pred, ylim = (0,6))
    display(plt)
    return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,p)->loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, p)

result_ode = Optimization.solve(optprob, PolyOpt(), callback = callback, maxiters = 100)

remade_solution = solve(remake(prob, p = result_ode.u), Tsit5(), saveat=tsteps)
plot(remade_solution, ylim = (0,6))