using Bridge
using Plots

# Define a diffusion process 

struct OrnsteinUhlenbeck  <: ContinuousTimeProcess{Float64}
    β::Float64 # drift parameter (also known as inverse relaxation time)
    σ::Float64 # diffusion parameter
    function OrnsteinUhlenbeck(β::Float64, σ::Float64)
        #isnan(β) || β > 0. || error("Parameter λ must be positive.")
        #isnan(σ) || σ > 0. || error("Parameter σ must be positive.")
        new(β, σ)
    end
end

#define drift and diffusion coefficient



#simulate ornstein uhlenbeck using Euler scheme

W = sample(0:0.01:10, Wiener())
X = solve(EulerMaruyama(), 0.1, W, OrnsteinUhlenbeck(0.0,1.0))
Plots.plot(X, level = "X")
Plots.plot!(W.tt, W.yy)
X


# This allows us to sample one brownian Bridge, we want multiple
t = 2
n = 100
dt = t/n

x0 = 1.

P = Wiener()

typeof(P)

W = sample(range(0., stop =t, length=n), P)
W = sample(0:dt:t,P)

println(W.tt)

Plots.plot(W.tt, W.yy)

X = solve(EulerMaruyama(), x0, W, P)
Plots.plot!(X.tt, X.yy)

#Sample a Brownian Bridge ending in v at time s
s = 1.
v = 0.

B = sample(0:dt:s, WienerBridge(s,v))
Plots.plot(W.tt, W.yy, color="blue")
Plots.plot!(B.tt, B.yy, color = "red")

#Make a loop that repeats n bridge samples, and samples them at specific time

N = 1:1000
Xhat = Vector{Float64}()
anc = 0.5

for n in N
    B = sample(0:dt:s, WienerBridge(s,v))
    plot!(B.tt, B.yy)
    idx = findall(x -> x == anc, B.tt)
    val = B.yy[idx][1]
    push!(Xhat, val)
end



using Statistics

#Kernel plots here could be useful

function Phylo_Bridge(start, fin, fin_time, anc, dt, samples)
    N = 1:samples
    Xhat = Vector{Float64}()

    for n in N
        B = sample(0:dt:fin_time, WienerBridge(fin_time,fin), start)
        idx = findall(x -> x == anc, B.tt)
        val = B.yy[idx][1]
        push!(Xhat, val)
    end

return Xhat
end

test = Phylo_Bridge(5., 4., 10., 3., 0.1, 100.)



import Bridge: b, σ, a, transitionprob
Bridge.b(t,x, P::OrnsteinUhlenbeck) = -P.β*x
Bridge.σ(t, x, P::OrnsteinUhlenbeck) = P.σ
Bridge.a(t, x, P::OrnsteinUhlenbeck) = P.σ^2

# simulate OrnsteinUhlenbeck using Euler scheme
W = sample(0:0.01:10, Wiener{Float64}()) 
X = solve(EulerMaruyama(), 0.1, W, OrnsteinUhlenbeck(20., 1.))

plot(X.tt, X.yy)

# define transition density
transitionprob(s, x, t, P::OrnsteinUhlenbeck) = Normal(x*exp(-P.β*(t-s)), sqrt((0.5P.σ^2/P.β) *(1-exp(-2*P.β*(t-s)))))

# plot likelihood of β 
LL = [(β, llikelihood(X, OrnsteinUhlenbeck(β, 1.))) for β in 1.:30.]
for (β, ll) in LL
    println("β $β loglikelihood ", ll )
end
plot(Float64[β for (β, ll) in LL], Float64[ll for (β, ll) in LL])

function Phylo_Bridge_Plot(start, final, dt, t, n)
    bridge_data = DataFrame()
    for i in 1:n
        startval = sample(start)
        finalval = sample(final)
        B = sample(0:dt:t, WienerBridge(t, finalval), startval)
        bridge_data[!, string(i)] = Vector(B.yy)
    end
    return(bridge_data)
end