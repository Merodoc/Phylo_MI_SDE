#Constants

D = 1e-2 #Diffusion Coefficient
dt = 1e-3 #Time Step
Np = 5e2 # No. Particles
NT = 1e3 # No. Time steps

x = zeros(Int(Np), Int(NT))
y = copy(x) 
z = copy(x)

x2 = copy(x)
y2 = copy(x)
z2 = copy(x)

using Random, Distributions
Random.seed!(123)

#=
 Populate Particle Movements
 sqrt(2Ddt)W(t) is the zero mass limit of LE
 =#

 for k in collect(2:Int(NT))
    # Langevin Equation
    x[:,k] = x[:,k-1] + sqrt(2D*dt)*randn(Int(Np),1)
    y[:,k] = y[:,k-1] + sqrt(2D*dt)*randn(Int(Np),1)
    z[:,k] = z[:,k-1] + sqrt(2D*dt)*randn(Int(Np),1)
    # Brownian Motion
    x2[:,k] = x2[:,k-1] + rand(Normal(0,0.1),Int(Np))
    y2[:,k] = y2[:,k-1] + rand(Normal(0,0.1),Int(Np))
    z2[:,k] = z2[:,k-1] + rand(Normal(0,0.1),Int(Np))
    # FPE
 end
