using Bridge
using StaticArrays

struct OU <: ContinuousTimeProcess{Float64}
    μ::Float64
end
Bridge.b(s, x, P::OU) = -P.μ*x
Bridge.σ(s, x, P::OU) = I

solve(EulerMaruyama(), 1.0, sample(0:0.1:10, Wiener()), OU(1.4))

SV = SVector{3,Float64}

tt = collect(0:0.01:60); length(tt)
U = SamplePath(tt, zeros(SV, length(tt)))

struct LogoP <: ContinuousTimeProcess{SV}
end
Bridge.b(s, x, ::LogoP) = b(s, x) # drift
Bridge.σ(s, x, ::LogoP) = 0.1I
Bridge.a(s, x, ::LogoP) = 0.1*0.1*I # diffusivity

X = copy(U)

z1, z2, z3 = zz = [SV(0,0,0), SV(-1,-1/sqrt(3),0), SV(-1, 1/sqrt(3), 0)]
u = z3

solve!(Bridge.EulerMaruyama(), X, u, sample(tt, Wiener{SV}()), LogoP() );