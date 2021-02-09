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

# Onto 2D stuff! Also, I'm not good at updating these comments, so sorry to anyone who reads these in the future (except you Dave, you deserve it for forgetting to update the comments.)

# This is what I think is typically called the "multi-grid method", but I understand it far more intuitively from this discretization viewpoint.

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
u = ones(ny, nx) # X-velocity
u[round(Int, 0.5/dy):round(Int, 1/dy+1), round(Int, 0.5/dx):round(Int, 1/dx+1)] .= 2

v = ones(ny, nx) # Y-velocity
v[round(Int, 0.5/dy):round(Int, 1/dy+1), round(Int, 0.5/dx):round(Int, 1/dx+1)] .= 2



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
    local vn = copy(v)
    for j in 2:ny-1
        for i in 2:nx-1
            u[i, j] = un[i, j] - un[i, j] * dt / dx * (un[i, j] - un[i-1, j]) - vn[i, j] * dt / dy * (un[i, j] - un[i, j-1]) # 2D, non-linear convection.
            v[i, j] = vn[i, j] - un[i, j] * dt / dx * (vn[i, j] - vn[i-1, j]) - vn[i, j] * dt / dy * (vn[i, j] - vn[i, j-1]) # Also 2D, non-linear convection.

            # Apply boundary conditions
            u[1, :] .= 1
            u[end, :] .= 1
            u[:, 1] .= 1
            u[:, end] .= 1
            v[1, :] .= 1
            v[end, :] .= 1
            v[:, 1] .= 1
            v[:, end] .= 1
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
