using Bridge
using Plots

# Define a diffusion process 

struct OrnsteinUhlenbeck <: ContinuousTimeProcess{Float64}
    β::Float64 #drift parameter 
    σ::Float64 #diffusion parameter
end

#define drift and diffusion coefficient
Bridge.b(t,x,P::OrnsteinUhlenbeck) = -P.β*x
Bridge.σ(t,x,P::OrnsteinUhlenbeck) = P.σ

#simulate ornstein uhlenbeck using Euler scheme

W = sample(0:0.01:10, Wiener())
X = solve(EulerMaruyama(), 0.1, W, OrnsteinUhlenbeck(2.0,1.0))
plot(X, level = "X")
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

plot(W.tt, W.yy)

X = solve(EulerMaruyama(), x0, W, P)

#Sample a Brownian Bridge ending in v at time s
s = 1.
v = 0.

B = sample(0:dt:s, WienerBridge(s,v))
plot(W.tt, W.yy, color="blue")
plot!(B.tt, B.yy, color = "red")

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