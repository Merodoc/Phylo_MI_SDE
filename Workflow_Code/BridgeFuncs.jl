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
plot!(X2)

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

B = sample(0:0.01:1, WienerBridge(5.,5.), 5.)
P = BridgeProp(OrnsteinUhlenbeck(1.,0.1), 0:0.01:1, (5., 5.), 1.)
P2 = BridgeProp(BMDrift(0.5, 1.), 0:0.01:10, (0.,0.),1.)
P3 = BridgeProp(OUMean(1.,5., 0.1), 0:0.01:1, (5.,5.),1.)
P4 = BridgeProp(CIR(0.5,0.,1.), 0:0.01:10,(0.,0.),1.)
P5 = BridgeProp(WF(0.5,0.,1.), 0:0.01:10,(0.,0.),1.)


X = solve(EulerMaruyama(), 0.1,B, P)
X2 = solve(EulerMaruyama(), 0.1, B, P2)
X3 = solve(EulerMaruyama(), B, P3)
X4 = solve(EulerMaruyama(), 0.1, B, P4)
X5 = solve(EulerMaruyama(), 0.1, B, P5)

plot(X)
plot(X2)
plot(X3)
plot!(X4)
plot!(X5)