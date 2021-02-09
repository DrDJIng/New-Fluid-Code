using LinearAlgebra
using Plots
# plotly()

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

# Onto Burger's equation! Combining both advection AND diffusion! Combining the two yields:
# u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + nu * dt / (dx * dx) * (un[i+1] - 2 * un[i] + un[i-1])

# Though we want to add some boundary conditions here, because we want to explore some interesting concepts with Burger's equation.
# Initla conditions are:
# u = -2nu/phi dphi/dx +4
# phi = exp(-x^2/4nu) + exp(-(x-2pi)^2/4nu)

# Has an analytic solution of:
# u = -2nu/phi dphi/dx +4
# phi = exp(-(x - 4t)^2/(4nu(t+1))) + exp(-(x-4t-2pi)^2/(4nu(t+1)))

# Boundary condition is: u(0) = u(2pi) -> Periodic!

nx = 101
ny = 101
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
nt = 100 # Number of timesteps
c = 1
sigma = 0.2 # Sigma is randomly defined later, then we're told to ignore it for now. K. Probably some kind of discretised stepping rate.
dt = sigma * dx # dt is the amount of time each step, i.e. delta t

t = 0 # Initial time
x = range(0, 2, length = nx)
y = range(0, 2, length = ny)

# Set up initial conditions. We want 1 everywhere, except where there are twos (0.5 ≤ (x,y) ≤ 1)
u = ones(ny, nx)
u[round(Int, 0.5/dy):round(Int, 1/dy+1), round(Int, 0.5/dx):round(Int, 1/dx+1)] .= 2



uini = copy(u)
p1 = contour(x, y, uini,
            xlim = [0, 2],
            ylim = [0, 2],
            fill = true,
            aspect_ratio = :equal,
            cbar = false)

# For every element in u, we need to perform the discretization
# anim = @animate for n in 1:nt
for n in 1:nt # This n actually isn't doing anything at the moment. I guess we *could* use it to store time information, but atm we just lose that.
    local un = copy(u)
    for j in 2:ny-1
        for i in 2:nx-1
            u[i, j] = un[i, j] - c * dt / dx * (un[i, j] - un[i-1, j]) - c * dt / dy * (un[i, j] - un[i, j-1]) # Linear convection. c -> Linear!

            # Apply boundary conditions
            u[1, :] .= 1
            u[end, :] .= 1
            u[:, 1] .= 1
            u[:, end] .= 1
        end
    end
    # plot(uini)
    # plot!(u)
end

p2 = contour(x, y, u,
            aspect_ratio = :equal,
            xlim = [0, 2],
            ylim = [0, 2],
            cbar = false)
# gif(anim, "sawtooth.gif", fps = 30)
plot(p1, p2, fill = true, layout = (1, 2))
