using DiffusionMCMC, ExtensibleMCMC, GuidedProposals
using DiffusionDefinition, ObservationSchemes
const DD = DiffusionDefinition
const eMCMC = ExtensibleMCMC

using StaticArrays, Plots
using OrderedCollections

@load_diffusion FitzHughNagumo

# generate some data
θ = [0.1, -0.8, 1.5, 0.0, 0.3]
P = FitzHughNagumo(θ...)
tt, y1 = 0.0:0.0001:10.0, @SVector [-1.0, -1.0]
X = rand(P, tt, y1)
data = map(
    x->(x[1], x[2][1] + 0.1randn()),
    collect(zip(X.t, X.x))[1:1000:end]
)[2:end]

# and let's see how they look
plot(X, Val(:vs_time), size=(800, 300))
scatter!(getindex.(data, 1), getindex.(data, 2), label="data")

recording = (
    P = P,
    obs = load_data(
        ObsScheme(
            LinearGsnObs(
                0.0, (@SVector [0.0]);
                L=(@SMatrix [1.0 0.0]), Σ=(@SMatrix [0.01])
            )
        ),
        data
    ),
    t0 = 0.0,
    x0_prior = KnownStartingPt(y1),
)
observs = AllObservations()
add_recording!(observs, recording)
init_obs, _ = initialize(observs)

# auxiliary law
@load_diffusion FitzHughNagumoAux

# main definition of an MCMC algorithm
mcmc_params = (
    mcmc = MCMC(
        [
            # sampling of paths with a precond. Crank–Nicolson parameter ρ=0.96
            PathImputation(0.96, FitzHughNagumoAux),
        ];
        backend=DiffusionMCMCBackend(),
    ),
    num_mcmc_steps = Integer(1e3),
    data = init_obs,
    θinit = nothing, # no inference
)

mcmc_kwargs = (
    path_buffer_size = 10,
    dt = 0.001,
)

# run the MCMC
glob_ws, loc_ws = run!(mcmc_params...; mcmc_kwargs...)

xs = glob_ws.XX_buffer.XX

p = plot(size=(800, 300))
for i in 1:10
    plot!(p, xs[i][1], Val(:vs_time), color=["red" "steelblue"], label=false, alpha=0.4)
end
scatter!(p, getindex.(data, 1), getindex.(data, 2), label="data")
display(p)