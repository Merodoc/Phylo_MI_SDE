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

θ = [0.0, 3.2, #==# 1.0 #==#]
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
DD.b(t, x, P::OrnsteinUhlenbeckAux) = DD.B(t,P)*x + DD.β(t,P)
DD.a(t, P::OrnsteinUhlenbeckAux) = DD.σ(t,P)*DD.σ(t,P)
DD.constdiff(P::OrnsteinUhlenbeckAux) = true

using ObservationSchemes, StaticArrays
t, v = 12.52, @SVector [3.06]
obs = LinearGsnObs(t, v; full_obs=true)

dt = 0.001

tt=0.0:dt:t
P = GuidProp(tt, P_target, OrnsteinUhlenbeckAux, obs)

x0 = @SVector[2.8]

X, W, Wnr = rand(P, x0)

plot(X, Val(:vs_time))

# Try Parameter Inference Hardcode

function customkernel(θ, s::Symbol, scale=0.1)
	θ° = deepcopy(θ)
	θ°[s] += 2.0*scale*(rand()-0.5)
	θ°
end

recording = (
    P = P_target,
	obs = load_data(
        ObsScheme(
            LinearGsnObs(t, v; full_obs=true)),
            [12.52], [v]),
	t0 = 0.0,
	x0_prior = KnownStartingPt(x0)
)



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

DD.const_parameter_names(::Type{<:OrnsteinUhlenbeck}) = (:θ, :μ)
DD.const_parameter_names(::Type{<:OrnsteinUhlenbeckAux}) = (:θ, :μ, :t0, :T, :vT)

paths, θθ = simple_inference(
	OrnsteinUhlenbeckAux, recording, 0.001, Dict(:σ=>1.0); ρ=0.5, num_steps=10^4
)

θθ
plot(θθ)

p = plot(size=(1400, 800))
for path in paths[end-10:end]
	for i in eachindex(path)
		plot!(p, path[i], Val(:vs_time), alpha=0.4, label="", color=["red" "steelblue"])
	end
end
display(p)



