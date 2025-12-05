using Plots, DiffusionDefinition, StaticArrays, GuidedProposals

const DD = DiffusionDefinition;

#Define Target law

@diffusion_process OrnsteinUhlenbeck begin
    :dimensions
    process --> 1
    wiener --> 1
    
    :parameters
    (θ, μ, σ) --> Float64
end

DD.b(t, x, P::OrnsteinUhlenbeck) = P.θ*(P.μ - x)
DD.σ(t, x, P::OrnsteinUhlenbeck) = P.σ

DD.default_type(::OrnsteinUhlenbeck) = Float64
DD.default_wiener_type(::OrnsteinUhlenbeck) = Float64

θ = [2.0, 3.2, 0.1]
P_target = OrnsteinUhlenbeck(θ...)

#Define Auxiliary law

@diffusion_process OrnsteinUhlenbeckAux begin
    :dimensions
    process --> 1
    wiener --> 1
    
    :parameters
    (θ, μ, σ) --> Float64

    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> SVector
end

DD.B(t, P::OrnsteinUhlenbeckAux) = @SMatrix[0.0]
DD.β(t, P::OrnsteinUhlenbeckAux) = @SVector[0.0]
DD.σ(t, P::OrnsteinUhlenbeckAux) = P.σ
constdiff(::OrnsteinUhlenbeckAux) = true
DD.b(t, x, P::OrnsteinUhlenbeckAux) = B(t,P)*x + β(t,P)
DD.a(t, P::OrnsteinUhlenbeckAux) = DD.σ(t,P)*DD.σ(t,P)


using ObservationSchemes, StaticArrays
t, v = 12.52, @SVector [3.06]
obs = LinearGsnObs(t, v; full_obs=true)

dt = 0.001

tt=0.0:dt:t
P = GuidProp(tt, P_target, OrnsteinUhlenbeckAux, obs)

@load_diffusion LotkaVolterraAux

x0 = @SVector[2.8]

X, W, Wnr = rand(P, x0)

plot(X, Val(:vs_time))