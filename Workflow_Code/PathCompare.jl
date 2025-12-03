include("SDEDefs.jl")

θ = [0.0, 3.2, #==# 1.0 #==#]
P_target = OrnsteinUhlenbeck(θ...)
P_target2 = CIR(θ...)

X = rand(P_target2, tt, 1.0)
using ObservationSchemes, StaticArrays
t, v = 12.52, @SVector [3.06]
obs = LinearGsnObs(t, v; full_obs=true)

dt = 0.001

tt=0.0:dt:t
P = GuidProp(tt, P_target, OrnsteinUhlenbeckAux, obs)
P2 = GuidProp(tt, P_target2, CIRAux, obs)

x0 = @SVector[2.8]

X, W, Wnr = rand(P, x0)
X2, W, Wnr = rand(P2,x0)

using Plots
plot(X, Val(:vs_time))
plot(X2, Val(:vs_time))

success, ll = GP.rand!(P, X2, W, Val(:ll), x0; Wnr=Wiener())


function simple_smoothing(P, y1)
	X, W, Wnr = rand(P, y1)
	X°, W° = trajectory(P)

	ll = loglikhd(P, X)
	paths = []

	for i in 1:10^4
		_, ll° = GP.rand!(P, X°, W°, Val(:ll), y1; Wnr=Wnr)
		if rand() < exp(ll°-ll)
			X, W, X°, W° = X°, W°, X, W
			ll = ll°
		end
		i % 400 == 0 && append!(paths, [deepcopy(X)])
	end
	paths
end
paths = simple_smoothing(P2, x0)

using Plots, Colors
cm = colormap("RdBu")
kwargs = (alpha=0.3, label="")
p = plot(paths[1], Val(:vs_time); color=cm[1], kwargs...)
for (i,x) in enumerate(paths[2:end])
	plot!(p, x, Val(:vs_time); color=cm[4*i], kwargs...)
end
display(p)
