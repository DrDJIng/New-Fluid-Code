using LinearAlgebra
using Plots
plotly()

# Basic equations
# Velocity field
# Dv/Dt = -1/ρ ∇p + ν∇^2v + F (i.e. Navier-Stokes)
# Density
# Dρ/Dt + ρ∇⋅v = 0
# Composition
# DS/Dt = S{dotted}

# Probably also need thermodynamic equations as well???  Let's just stick to velocity/density for the moment.

# So, the change in velocity is given by the rate of change of the pressure, and the jerk of the velocity field.

# Basic algorithm should go something like:
# Set up initial matrices
# Iteratively solve integrals across matrices
# Update initial matrices
# Repeat for more timesteps

# First, let's build up from 1D, and go from there.

nx = 61
dx = 2 / (nx - 1)
nt = 2 # Number of timesteps
dt = 0.025 # dt is the amount of time each step, i.e. delta t
c = 1 # Wavespeed

# Set up initial conditions. We want 1 everywhere, except where there are twos (0.5 ≤ x ≤ 1)
u = ones(1, nx)
u[round(Int, 0.5/dx):round(Int, 1/dx+1)] .= 2

plot(u')

# For every element in u, we need to perform the discretization
for n in 1:nt
    local u_n = copy(u)
    for i in 2:nx
        u[i] = u_n[i] - u_n[i] * dt / dx * (u_n[i] - u_n[i-1]) # Non-linear term added here (c -> u_n[i]).
    end
end

plot!(u')
