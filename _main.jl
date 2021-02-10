using LinearAlgebra
using Plots

# Basic equations
# Velocity field
# Dv/Dt = -1/ρ ∇p + ν∇^2v + F (i.e. Navier-Stokes)
# Density
# Dρ/Dt + ρ∇⋅v = 0
# Composition
# DS/Dt = S{dotted}

# Probably also need thermodynamic equations as well???  Let's just stick to velocity/density for the moment.

# So, the change in velocity is given by the rate of change of the p, and the jerk of the velocity field.

# Basic algorithm should go something like:
# Set up initial matrices
# Iteratively solve integrals across matrices
# Update initial matrices
# Repeat for more timesteps

# Onto 2D stuff! Also, I'm not good at updating these comments, so sorry to anyone who reads these in the future (except you Dave, you deserve it for forgetting to update the comments.)

# Poisson equation is just the Laplace equation, but with a source term.

# nx = 41
# ny = 41
# dx = 2 / (nx - 1)
# dy = 2 / (ny - 1)
# nt = 120 # Number of timesteps
# nu = 0.01
# sigma = 0.0009 # Sigma is randomly defined, then we're told to ignore it for now. K. Probably some kind of discretised stepping rate.
# dt = sigma * dx * dy / nu # dt is the amount of time each step, i.e. delta t

# t = 0 # Initial time
# x = range(0, 2, length = nx)
# y = range(0, 2, length = ny)

# Set up initial conditions. We want 1 everywhere, except where there are twos (0.5 ≤ (x,y) ≤ 1)
# u = ones(ny, nx) # X-velocity
# u[round(Int, 0.5/dy):round(Int, 1/dy+1), round(Int, 0.5/dx):round(Int, 1/dx+1)] .= 2
#
# v = ones(ny, nx) # Y-velocity
# v[round(Int, 0.5/dy):round(Int, 1/dy+1), round(Int, 0.5/dx):round(Int, 1/dx+1)] .= 2
#
# uini = copy(u)
# p1 = contour(x, y, uini,
#             xlim = [0, 2],
#             ylim = [0, 2],
#             fill = true,
#             clim = (1, 2),
#             aspect_ratio = :equal,
#             #cbar = false
#             )


# ufin = diffuse!(nt, u)

# This is back in as explicit code, because fuck functions, I guess?

# for n in 1:nt # This n actually isn't doing anything at the moment. I guess we *could* use it to store time information, but atm we just lose that.
#     local un = copy(u)
#     local vn = copy(v)
#     for j in 2:ny-1
#         for i in 2:nx-1
#             # u[i, j] = un[i, j] + nu * dt / (dx * dx) * (un[i+1, j] - 2 * un[i, j] + un[i-1, j]) + nu * dt / (dy * dy) * (un[i, j+1] - 2 * un[i, j] + un[i, j-1]) # 2D diffusion
#             # u[i, j] = un[i, j] - un[i, j] * dt / dx * (un[i, j] - un[i-1, j]) - vn[i, j] * dt / dy * (un[i, j] - un[i, j-1]) # 2D, non-linear convection.
#             # v[i, j] = vn[i, j] - un[i, j] * dt / dx * (vn[i, j] - vn[i-1, j]) - vn[i, j] * dt / dy * (vn[i, j] - vn[i, j-1]) # Also 2D, non-linear convection.
#             u[i, j] = (un[i, j]
#                     - dt / dx * un[i, j] * (un[i, j] - un[i-1, j])
#                     - dt / dy * vn[i, j] * (un[i, j] - un[i, j-1])
#                     + nu * dt / (dx * dx) * (un[i+1, j] - 2 * un[i, j] + un[i-1, j])
#                     + nu * dt / (dy * dy) * (un[i, j+1] - 2 * un[i, j] + un[i, j-1])
#                     )
#
#             v[i, j] = (vn[i, j]
#                     - dt / dx * un[i, j] * (vn[i, j] - vn[i-1, j])
#                     - dt / dy * vn[i, j] * (vn[i, j] - vn[i, j-1])
#                     + nu * dt / (dx * dx) * (vn[i+1, j] - 2 * vn[i, j] + vn[i-1, j])
#                     + nu * dt / (dy * dy) * (vn[i, j+1] - 2 * vn[i, j] + vn[i, j-1])
#                     )
#
#             # Apply boundary conditions
#             u[1, :] .= 1
#             u[end, :] .= 1
#             u[:, 1] .= 1
#             u[:, end] .= 1
#             v[1, :] .= 1
#             v[end, :] .= 1
#             v[:, 1] .= 1
#             v[:, end] .= 1
#         end
#     end
#     # plot(uini)
#     # plot!(u)
# end


# p2 = contour(x, y, u,
#             aspect_ratio = :equal,
#             xlim = [0, 2],
#             ylim = [0, 2],
#             clim = (1, 2),
#             #cbar = false
#             )
# # gif(anim, "sawtooth.gif", fps = 30)
# plot(p1, p2, fill = true, layout = (1, 2))

function laplace2D!(p, b, y, nx::Int, ny::Int, dx::Float64, dy::Float64, err::Float64)
        current_err = 1

        while current_err > err
                local pn = copy(p)
                for j = 2:ny - 1
                        for i = 2:nx - 1
                                p[i, j] = ((dy * dy * (pn[i+1, j] + pn[i-1, j]) + dx * dx * (pn[i, j+1] + pn[i, j-1]) - b[i, j] * dx * dx * dy * dy)
                                        / (2 * (dx * dx + dy * dy))
                                        ) # Laplace is deceptively simple, isn'tit?
                        end
                end
                p[:, 1] .= 0 # p = 0 @ x = 0
                p[:, nx-1] .= 0 # p = 0 @ x = 2
                p[1, :] .= 0 # dp/dy = 0 @ y = 0
                p[ny-1, :] .= 0 # dp/dy = 0 @ y = 1

                current_err = (sum(abs.(p) - abs.(pn))) / sum(abs.(pn))
        end
        return nothing
end

nx = 31
ny = 31
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
x = range(0, 2, length = nx)
y = range(0, 1, length = ny)
p = zeros(ny, nx)

b = zeros(size(p))
b[floor(Int, ny/4), floor(Int, nx/4)] = 100
b[floor(Int, 3*ny/4), floor(Int, 3*nx/4)] = -100

# p1 = contour(x, y, p,
#             # aspect_ratio = :equal,
#             xlim = [0, 2],
#             ylim = [0, 1],
#             # clim = (1, 2),
#             #cbar = false
#             )

laplace2D!(p, b, y, nx, ny, dx, dy, 1e-4)

p2 = contour(x, y, p,
            # aspect_ratio = :equal,
            xlim = [0, 2],
            ylim = [0, 1],
            # clim = (-0.25, 0.25),
            #cbar = false
            )

plot(p2, fill = true)
