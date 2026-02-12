using Bridge

struct OrnsteinUhlenbeck <: ContinuousTimeProcess{Float64}
    β::Float64 #drift parameter 
    σ::Float64 #diffusion parameter
end

Bridge.b(t,x,P::OrnsteinUhlenbeck) = -P.β*x
Bridge.σ(t,x,P::OrnsteinUhlenbeck) = P.σ

using Random
Random.seed!(1)
W = sample(0:0.1:1, Wiener())

X = Bridge.solve(Euler(), 0.1, W, OrnsteinUhlenbeck(20.0, 1.0));

X

struct BMDrift <: ContinuousTimeProcess{Float64}
    β::Float64 #drift
    σ::Float64 #Diffusion
end

Bridge.b(t, x, P::BMDrift) = P.β
Bridge.σ(t, x, P::BMDrift) = P.σ

struct OUMean <: ContinuousTimeProcess{Float64}
    α::Float64 #drift
    μ::Float64 #Mean
    σ::Float64 #Diffusion
end

Bridge.b(t, x, P::OUMean) = P.α*(P.μ-x)
Bridge.σ(t,x,P::OUMean) = P.σ

struct CIR <: ContinuousTimeProcess{Float64}
    α::Float64 #drift
    μ::Float64 #Mean
    ϵ::Float64 #Diffusion
end

Bridge.b(t, x, P::CIR) = P.α*(P.μ-x)
Bridge.σ(t,x,P::CIR) = (P.ϵ*x)^1/2

struct WF <: ContinuousTimeProcess{Float64}
    α::Float64 #drift
    μ::Float64 #Mean
    σ::Float64 #Diffusion
end

Bridge.b(t, x, P::WF) = P.α*(P.μ-x)
Bridge.σ(t,x,P::WF) = (P.σ*x*(1-x))^1/2





X5 = Bridge.solve(Euler(), 0.1, W, CIR(20.0, 0.0, 1.0));

X   

using Plots

plot(X)

function Phylo_BridgeSDEtest(start, fin, fin_time, anc, dt, samples, SDE = OrnsteinUhlenbeck(0.0, 1.0))
    #start = start value
    #fin = final value
    #fin_time = total time
    #anc = time point of the ancestor node
    #samples = number of repeats
    N = 1:samples
    Xhat = Vector{Float64}()


    # Define the OU Process


    #n is no.times we are running the bridge... I can change this (run more trees)

    for n in N
        # Solve an OU bridge beteen start and fin
        B = sample(0:dt:fin_time, WienerBridge(fin_time,fin), start)
        #
        X = solve(EulerMaruyama(), 0.1, B, SDE)
        idx = findall(x -> x == anc, X.tt)
        val = X.yy[idx][1]
        push!(Xhat, val)
    end

return Xhat
end




B = sample(0:0.01:5, WienerBridge(5.,5.), 5.)
P2 = BridgeProp(BMDrift(0.5, 0.1), 0:0.01:5, (5.,5.),5.)
P3 = BridgeProp(OUMean(1.,5., 0.1), 0:0.01:5, (5.,5.),5.)
P4 = BridgeProp(CIR(1,5.,0.1), 0:0.01:5,(5.,5.),5.)
P5 = BridgeProp(WF(1,5.,.1), 0:0.01:5,(5.,5.),5.)


X2 = solve(EulerMaruyama(), B, P2)
X3 = solve(EulerMaruyama(), B, P3)
X4 = solve(EulerMaruyama(), B, P4)
X5 = solve(EulerMaruyama(), B, P5)

plot(X2)
plot!(X3)
plot!(X4)
plot!(X5)


p = plot()
OUspread = Float64[]

using Statistics
using StatsKit
using KernelDensity

Models = [P3, P4, P5]
p = plot()
for model in Models
    midpoints = Float64[]
    for i in 1:1000
        B = sample(0:0.01:5, WienerBridge(5.,5.), 5.)
        Y = solve(EulerMaruyama(), 5., B, model)
        push!(midpoints, Y.yy[250])
    end
    
    kernel = kde(midpoints)
    plot!(p, kernel.x, kernel.density)
end

display(p)




kernel = kde(OUspread)
plot(kernel.x, kernel.density)
Y.yy[250]

