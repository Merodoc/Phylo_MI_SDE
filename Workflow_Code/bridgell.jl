using Bridge, Distributions, StaticArrays
using Plots
using LinearAlgebra

struct OrnsteinUhlenbeck  <: ContinuousTimeProcess{Float64}
    β::Float64 # drift parameter (also known as inverse relaxation time)
    σ::Float64 # diffusion parameter
    function OrnsteinUhlenbeck(β::Float64, σ::Float64)
        #isnan(β) || β > 0. || error("Parameter λ must be positive.")
        #isnan(σ) || σ > 0. || error("Parameter σ must be positive.")
        new(β, σ)
    end
end

# define drift and sigma of OrnsteinUhlenbeck
import Bridge: b, σ, a, transitionprob
Bridge.b(t,x, P::OrnsteinUhlenbeck) = -P.β*x
Bridge.σ(t, x, P::OrnsteinUhlenbeck) = P.σ
Bridge.a(t, x, P::OrnsteinUhlenbeck) = P.σ^2

# simulate OrnsteinUhlenbeck using Euler scheme
W = sample(0:0.01:12.52, Wiener{Float64}())

dt = 0.01
s = 12.52
v = 3.06
x0 = 2.81
B = sample(0:dt:s, WienerBridge(s,v), x0)
X = solve(EulerMaruyama(), 0.1, B, OrnsteinUhlenbeck(10., 1.))

plot(X.tt, X.yy)

# define transition density
transitionprob(s, x, t, P::OrnsteinUhlenbeck) = Normal(x*exp(-P.β*(t-s)), sqrt((0.5P.σ^2/P.β) *(1-exp(-2*P.β*(t-s)))))

# plot likelihood of β 
LL = [(β, llikelihood(X, OrnsteinUhlenbeck(β, 1.))) for β in 1.:30.]
for (β, ll) in LL
    println("β $β loglikelihood ", ll )
end
plot(Float64[β for (β, ll) in LL], Float64[ll for (β, ll) in LL])

LL = Vector{}
for β in 1.:30.
    X = solve(EulerMaruyama(), 0.1, B, OrnsteinUhlenbeck(β, 1.))
    ll = llikelihood(X, OrnsteinUhlenbeck(β, 1.))
    push!{LL, (β, ll)}
end