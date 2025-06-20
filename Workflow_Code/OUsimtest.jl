include("OUOpt.jl")

θ = [0.0, 3.2, #==# 1.0 #==#]
P_target = OrnsteinUhlenbeck(θ...)

using ObservationSchemes, StaticArrays
t, v = 12.52, @SVector [3.06]
obs = LinearGsnObs(t, v; full_obs=true)

dt = 0.001

tt=0.0:dt:t
P = GuidProp(tt, P_target, OrnsteinUhlenbeckAux, obs)

x0 = @SVector[2.8]

X, W, Wnr = rand(P, x0)

plot(X, Val(:vs_time))

recording = (
    P = P_target,
	obs = load_data(
        ObsScheme(
            LinearGsnObs(t, v; full_obs=true)),
            [12.52], [v]),
	t0 = 0.0,
	x0_prior = KnownStartingPt(x0)
)



DD.const_parameter_names(::Type{<:OrnsteinUhlenbeck}) = (:θ, :μ)
DD.const_parameter_names(::Type{<:OrnsteinUhlenbeckAux}) = (:θ, :μ, :t0, :T, :vT)

paths, θθ = simple_inference(
	OrnsteinUhlenbeckAux, recording, 0.001, Dict(:σ=>0.0); ρ=0.8, num_steps=10^4, ϵ = 0.1
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

using Statistics
using StatsKit
test = kde(θθ)

plot(test.x, test.density)

mean(θθ)

#Lets try and hit a known parameter

tt, y1 = 0.0:0.0001:10.0, @SVector [1.0]
X = rand(P_target, tt, y1)
data = map(
    x ->(x[1], x[2][1] + 0.1randn()),
    collect(zip(X.t, X.x))[1:1000:end]
)[2:end]

plot(X, Val(:vs_time), size=(800, 300))
scatter!(map(x->x[1], data), map(x->x[2], data), label="data")

recording = (
	P = P_target,
	obs = load_data(
		ObsScheme(
			LinearGsnObs(
				0.0, (@SVector [0.0]);
				L=(@SMatrix [1.0]), Σ=(@SMatrix [0.01])
			)
		),
		data
	),
	t0 = 0.0,
	x0_prior = KnownStartingPt(y1),
)


DD.const_parameter_names(::Type{<:OrnsteinUhlenbeck}) = (:μ)
DD.const_parameter_names(::Type{<:OrnsteinUhlenbeckAux}) = (:θ, :μ, :t0, :T, :vT)

paths, θθ = simple_inference(
	OrnsteinUhlenbeckAux, recording, 0.001, Dict(:σ=>1.0, :θ => 0.0); ρ=0.8, num_steps=10^4, ϵ = 0.3
)

θθ

θθ
mean(θθ)
plot(θθ)

p = plot(size=(1400, 800))
for path in paths[end-10:end]
	for i in eachindex(path)
		plot!(p, path[i], Val(:vs_time), alpha=0.4, label="", color=["red" "steelblue"])
	end
end
scatter!(map(x->x[1], data), map(x->x[2], data), label="data")
display(p)
