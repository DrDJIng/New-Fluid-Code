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

# Discretising the diffusion involves a second-order derivative. Taylor-expand with foward and backward difference
# schemes, then combine to cancel terms. Neglect O(x^4) terms.

# When done properly, gives d2u/dx2 = (ui+1 - 2ui + ui-1)/Delx2 + O(Delx^2)

nx = 41
dx = 2 / (nx - 1)
nt = 20 # Number of timesteps
nu = 0.01 # Viscosity
sigma = 0.2 # Sigma is randomly defined later, then we're told to ignore it for now. K. Probably some kind of discretised stepping rate.
dt = sigma * dx * dx / ν # dt is the amount of time each step, i.e. delta t

# Set up initial conditions. We want 1 everywhere, except where there are twos (0.5 ≤ x ≤ 1)
u = ones(1, nx)
u[round(Int, 0.5/dx):round(Int, 1/dx+1)] .= 2

plot(u')

# For every element in u, we need to perform the discretization
for n in 1:nt
    local u_n = copy(u)
    for i in 2:nx-1
        u[i] = u_n[i] + nu * dt / (dx * dx) * (u_n[i+1] - 2 * u_n[i] + u_n[i-1]) # Discretised diffusion added here. It's linear, but based on viscosity.
    end
end

plot!(u')
