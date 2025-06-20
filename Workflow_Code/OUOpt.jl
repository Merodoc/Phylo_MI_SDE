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



# Try Parameter Inference Hardcode

function customkernel(θ, s::Symbol, scale=0.1)
	θ° = deepcopy(θ)
	θ°[s] += 2.0*scale*(rand()-0.5)
	θ°
end




# Perform inference on a single parameter for the data in the `recording`, using
# Guided Proposals with the auxiliary law `AuxLaw`.
function simple_inference(AuxLaw, recording, dt, θ; ρ=0.5, num_steps=10^4, ϵ = 0.3)
    # -------------------------------------------------------------------------#
    #                          Initializations                                 #
    # -------------------------------------------------------------------------#
    # time-grids for the forward-simulation of trajectories                    #
    tts = OBS.setup_time_grids(recording, dt)                                  #
    # laws of guided proposals                                                 #
    PP = build_guid_prop(AuxLaw, recording, tts)                               #
    # laws of guided proposals for parameter proposals                         #
    PP° = deepcopy(PP)                                                         #
                                                                               #
    # starting point                                                           #
    # NOTE `rand` for `KnownStartingPt` simply returns the starting position   #
    y1 = rand(recording.x0_prior)                                              #
    # initialize the `accepted` trajectory                                     #
    XX, WW, Wnr = rand(PP, y1)                                                 #
    # initialize the containers for the `proposal` trajectory                  #
    XX°, WW° = trajectory(PP)                                                  #
                                                                               #
    ll = loglikhd(PP, XX)                                                      #
    paths = []                                                                 #
    imp_num_accpt = 0                                                          #
    param_num_accpt = 0                                                        #
    pname = collect(keys(θ))[1]                                                #
    θθ = Float64[θ[pname],]                                                    #
    # -------------------------------------------------------------------------#

    # MCMC
    for i in 1:num_steps
        # impute a path
        _, ll° = DD.rand!(PP, XX°, WW°, WW, ρ, Val(:ll), y1; Wnr=Wnr)

        # Metropolis–Hastings accept/reject step
        if rand() < exp(ll°-ll)
            XX, WW, XX°, WW° = XX°, WW°, XX, WW
            ll = ll°
            imp_num_accpt += 1
        end

        # update parameter s
        θ° = customkernel(θ, pname, ϵ)
        DD.set_parameters!(PP°, θ°)
        recompute_guiding_term!(PP°)
        _, ll° = GP.solve_and_ll!(XX°, WW, PP°, y1)

        if rand() < exp(ll°-ll) # uniform updates have no contribution to ll
            XX, PP, θ, XX°, PP°, θ° = XX°, PP°, θ°, XX, PP, θ
        ll = ll°
        param_num_accpt += 1
    end
    append!(θθ, [θ[pname]])

    # progress message
    if i % 100 == 0
        println(
            "$i. ll=$ll, $pname=$(θ[pname]), imp accpt rate: ",
            "$(imp_num_accpt/100), updt accpt rate: $(param_num_accpt/100)"
        )
        imp_num_accpt = param_num_accpt = 0
    end

    # save intermediate path for plotting
    i % 400 == 0 && append!(paths, [deepcopy(XX)])
    end
    paths, θθ
end

