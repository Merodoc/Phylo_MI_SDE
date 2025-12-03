using Plots, DiffusionDefinition, StaticArrays, GuidedProposals, ObservationSchemes

const DD = DiffusionDefinition;
const OBS = ObservationSchemes
const GP = GuidedProposals

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
DD.constdiff(P::OrnsteinUhlenbeck) = true
DD.default_type(::OrnsteinUhlenbeck) = Float64
DD.default_wiener_type(::OrnsteinUhlenbeck) = Float64



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
DD.b(t, x, P::OrnsteinUhlenbeckAux) = DD.B(t,P)*x + DD.β(t,P)
DD.a(t, P::OrnsteinUhlenbeckAux) = DD.σ(t,P)*DD.σ(t,P)
DD.constdiff(P::OrnsteinUhlenbeckAux) = true



@diffusion_process CIR begin
    :dimensions
    process --> 1
    wiener --> 1
    
    :parameters
    (θ, μ, σ) --> Float64
end

DD.b(t, x, P::CIR) = P.θ*(P.μ - x)
DD.σ(t, x, P::CIR) = sqrt(P.σ*x[1])
DD.constdiff(P::CIR) = true
DD.default_type(::CIR) = Float64
DD.default_wiener_type(::CIR) = Float64

@diffusion_process CIRAux begin
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

DD.B(t, P::CIRAux) = @SMatrix[0.0]
DD.β(t, P::CIRAux) = @SVector[0.0]
DD.σ(t, P::CIRAux) = P.σ
DD.b(t, x, P::CIRAux) = DD.B(t,P)*x + DD.β(t,P)
DD.a(t, P::CIRAux) = DD.σ(t,P)*DD.σ(t,P)
DD.constdiff(P::CIRAux) = true


@diffusion_process WF begin
    :dimensions
    process --> 1
    wiener --> 1
    
    :parameters
    (θ, μ, σ) --> Float64
end

DD.b(t, x, P::WF) = P.θ*(P.μ - x)
DD.σ(t, x, P::WF) = sqrt(P.σ*x[1]*(1-x[1]))
DD.constdiff(P::WF) = true
DD.default_type(::WF) = Float64
DD.default_wiener_type(::WF) = Float64

@diffusion_process WFAux begin
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

DD.B(t, P::WFAux) = @SMatrix[0.0]
DD.β(t, P::WFAux) = @SVector[0.0]
DD.σ(t, P::WFAux) = P.σ
DD.b(t, x, P::WFAux) = DD.B(t,P)*x + DD.β(t,P)
DD.a(t, P::WFAux) = DD.σ(t,P)*DD.σ(t,P)
DD.constdiff(P::WFAux) = true