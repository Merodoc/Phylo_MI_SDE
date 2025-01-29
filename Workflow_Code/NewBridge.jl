using Bridge

# Define diffusion Process
struct OrnsteinUhlenbeck <: ContinuousTimeProcess{Float64}
    β::Float64
    σ::Float64
    function OrnsteinUhlenbeck(β::Float64, σ::Float64)
        isnan(β) || β > 0. || error("Parameter β must be positive.")
        isnan(σ) || σ > 0. || error("Parameter σ must be positive.")
        new(β, σ)
    end
end

#Define drift and diffusion coefficient

#dependent drift
Bridge.b(t, x, P::OrnsteinUhlenbeck) = -P.β * x
#dispersion coefficient
Bridge.σ(t, x, P::OrnsteinUhlenbeck) = P.σ
#diffusion coefficient
Bridge.a(t, x, P::OrnsteinUhlenbeck) = P.σ^2

#Simulate OU 
#Generate Driving Brownian Motion W of the SDE
using Random
Random.seed!(1)
W = sample(0:0.1:1, WienerBridge(1.0, 0.0))

#Solve

X = Bridge.solve(Euler(), 0.1, W, OrnsteinUhlenbeck(20.0, 1.0));
