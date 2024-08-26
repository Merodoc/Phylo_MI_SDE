library(pracma)

# Constants

D = 1e-2  #Diffusion Coefficient
dt = 1e-3 #time step
Np = 5e2  #Number of Particles
NT = 1e3  #Number of time steps

x = matrix(0, Np, NT)
y = matrix(0, Np,NT)
z = matrix(0, Np,NT)

x2 = matrix(0, Np, NT)
y2 = matrix(0, Np,NT)
z2 = matrix(0, Np,NT)

# Populate Particle Movements
# sqrt(2*D*dt)*W(t) is the zero mass limit of LE
for (k in 2:NT) {
  x[,k] = x[,k-1] + sqrt(2*D*dt)*randn(Np,1)
  y[,k] = y[,k-1] + sqrt(2*D*dt)*randn(Np,1)
  z[,k] = z[,k-1] + sqrt(2*D*dt)*randn(Np,1)
}

# As this is zero mass limit, the particles can sit inside each other
# Once we add spheres we should be able to prevent them from hitting each other
# Compare to Brownian Motion

for (k in 2:NT) {
  x2[,k] = x2[,k-1] + rnorm(Np, sd = 0.01)
  y2[,k] = y2[,k-1] + rnorm(Np, sd = 0.01)
  z2[,k] = z2[,k-1] + rnorm(Np, sd = 0.01)
}