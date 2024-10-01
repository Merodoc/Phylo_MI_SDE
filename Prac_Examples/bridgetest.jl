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
plot(X, lebel = "X")
X