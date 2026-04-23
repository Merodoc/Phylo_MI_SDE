using Plots, DiffusionDefinition, StaticArrays

const DD = DiffusionDefinition;

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

#sampling trajectory

tt, x0 = 0.0:0.01:12.52, 2.8
P = OrnsteinUhlenbeck(2.0, 3.2, 0.1)

X = rand(P, tt, x0)
#T = 12.52
#xT = 3.06
#x0 = 2.8

using Plots

gr()
plot(X, Val(:vs_time))


xT, T = 3.06, 12.52
B = trajectory(tt, Float64)
W = rand(Wiener(), tt, 0.0)
B.x .= x0 .+ W.x .+ tt./T.*(xT .- x0 .- W.x[end])

plot_kwargs = (color="steelblue", label="")
p = plot(B; plot_kwargs...)

for i in 1:30
    W = rand(Wiener(), tt, 0.0)
    B.x .= x0 .+ W.x .+ tt./T.*(xT .- x0 .- W.x[end])
    plot!(p, B; alpha=0.2, plot_kwargs...)
end
display(p)