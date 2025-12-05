using Bridge
using StaticArrays
import Bridge: b, σ, B, β, a, constdiff
using BridgeSDEInference

const ℝ = SVector{N,T} where {N,T}

"""
    RadialOU{T} <: ContinuousTimeProcess{ℝ{1,T}}

Struct defining Radial Ornstein Uhlenbeck processs
"""
struct RadialOU{T} <: ContinuousTimeProcess{ℝ{1,T}}
    η::T
    σ::T
    RadialOU(η::T, σ::T) where T = new{T}(η, σ)
end


b(t, x, P::RadialOU) = ℝ{1}(-P.η*x[1] + 0.5*P.σ^2/x[1])
σ(t, x, P::RadialOU) = ℝ{1}(P.σ^2)

domain(::RadialOU{T}) where T = LowerBoundedDomain((zero(T),), (1,))
constdiff(::RadialOU) = true
clone(P::RadialOU, θ) = RadialOU(θ...)
params(P::RadialOU) = [P.η, P.σ]


"""
    RadialOUAux{S} <: ContinuousTimeProcess{ℝ{1,S}}
"""
struct RadialOUAux{R, S1, S2} <: ContinuousTimeProcess{ℝ{1,R}}
    #trgtDomain::LowerBoundedDomain{S,1}
    η::R
    σ::R
    t::Float64
    u::S1
    T::Float64
    v::S2

    function RadialOUAux(η::R, σ::R, t::Float64, u::S1, T::Float64, v::S2) where {R, S1, S2}
        new{R, S1, S2}(η, σ, t, u, T, v)
    end

end

B(t, P::RadialOUAux) = @SMatrix[0.0]
β(t, P::RadialOUAux) = ℝ{1}(0.0)
σ(t, P::RadialOUAux) = ℝ{1}(P.σ^2)
dependsOnParams(::RadialOUAux) = (2,)
constdiff(::RadialOUAux) = true
b(t, x, P::RadialOUAux) = B(t,P)*x + β(t,P)
a(t, P::RadialOUAux) = σ(t,P)*σ(t,P)'
clone(P::RadialOUAux, θ) = RadialOUAux(P.trgtDomain, θ..., P.t, P.u, P.T, P.v)
params(P::RadialOUAux) = [P.η, P.σ]

using Random
Random.seed!(4)
θˣ = [1.5, 1.0] #  (α, β, γ, δ, σ1, σ2)
Pˣ = RadialOU(θˣ...)

x0 = ℝ{1}(2.8) # starting point
w0 = ℝ{1}(0.0) # starting point of driving Wiener process

dt, T =  1/5000, 12.52 # time grid
tt = 0.0:dt:T

dt = 0.01
s = 12.52
v = 3.06
xt = ℝ{1}(3.06)
x0 = ℝ{1}(2.8)
W = sample(0:dt:s, WienerBridge(s,v), 2.81)
X = solve(EulerMaruyama(), x0, W, Pˣ)

obs = W.yy
obs_times = [0.0, 12.52]
obs = []
obs = [x0, xt]
obs = [2.8, 3.06]
Σdiagel = 1.0
Σ = @SMatrix[Σdiagel]
L = @SMatrix[1.0 0.0]

θ_init = copy(θˣ)
Pˣ = RadialOU(θ_init...)

P̃ = [RadialOUAux(1.5, 1.0, 0.0, x0, 12.52, xt)]

model_setup = DiffusionSetup(Pˣ, P̃, PartObs())

set_observations!(model_setup, [L for _ in P̃], [Σ for _ in P̃], obs, obs_times)
set_imputation_grid!(model_setup, dt)
set_x0_prior!(model_setup, KnownStartingPt(x0))

initialise!(eltype(x0), model_setup, Vern7, false, NoChangePt)

obs_times