using LinearAlgebra
using AbstractPlotting
using GLMakie
AbstractPlotting.inline!(false)

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

# Finally getting into the nitty-gritty stuff. Proper Navier-Stokes! In a cavity, so just a closed box, but still on my way to actually ~understanding~
# the finite-difference method.

# I think I want two functions. A poisson pressure function, and a diffusion / advection function for the velocities

# Diffuse/advect should take in u/v, v/u, dx, dy, dt, rho, nu
# Should generalise so that I can call it twice with the u/v's swapped, since they're identical apart from that.
# Realised I can't do it the initial way I thought, because then the "guesses" wouldn't be the right guess for each run-through
# Well, it would be right for the first, then wrong for the second, as the first would already be iterated through time, and we need to consider
# guesses at the same time-step.

function moveFluid!(u, v, p, nx::Int, ny::Int, dx::Float64, dy::Float64, dt::Float64, rho::Float64, nu::Float64)
        local un = copy(u) # Make temporary copy of previous time step
        local vn = copy(v)
        for j = 2:ny - 1
                for i = 2:nx - 1
                        u[i, j] = (un[i, j] # Discretised Du/Dt = -1/ρ ∇p + ν∇^2u
                                - un[i, j] * (dt / dx) * (un[i, j] - un[i-1, j])
                                - vn[i, j] * (dt / dy) * (un[i, j] - un[i, j-1])
                                - dt / (rho * 2 * dx) * (p[i+1, j] - p[i-1, j])
                                + nu * (
                                        dt / (dx^2) * (un[i+1, j] - 2 * un[i, j] + un[i-1, j])
                                        + dt / (dy^2) * (un[i, j+1] - 2 * un[i, j] + un[i, j-1])
                                        )
                                )

                        v[i, j] = (vn[i, j] # Discretised Dv/Dt = -1/ρ ∇p + ν∇^2v
                                - un[i, j] * (dt / dx) * (vn[i, j] - vn[i-1, j])
                                - vn[i, j] * (dt / dy) * (vn[i, j] - vn[i, j-1])
                                - dt / (rho * 2 * dy) * (p[i, j+1] - p[i, j-1])
                                + nu * (
                                        dt / (dx^2) * (vn[i+1, j] - 2 * vn[i, j] + vn[i-1, j])
                                        + dt / (dy^2) * (vn[i, j+1] - 2 * vn[i, j] + vn[i, j-1])
                                        )
                                )

                                # Probably actually can put the un / vn in as an input. Yeah, that probably would have been the smarter thing to do
                                # Especially when this gets to 3D.
                end
        end

        # Apply boundary conditions.
        u[1, :] .= 0
        u[:, 1] .= 0
        u[end, :] .= 0
        u[:, end] .= 1 # Set horizontal movement on "the lid"
        v[1, :] .= 0
        v[end, :] .= 0
        v[:, 1] .= 0
        v[:, end] .= 0

        return nothing # Libbum tells me that this is faster than actually not returning anything.
end

# Poisson should take in p, u, v, dx, dy, dt. The version in the tutorial has a set number of iterations (nit), so maybe the err doesn't converge nicely?
# This is also where I will eventually put Jesse's nice stencils. Also just realised that iterative to equilibrium removes the divergence! Of course it does!
# I can be really thick sometimes.

function removeDivergenceFromPressure!(p, u, v, dx::Float64, dy::Float64, dt::Float64, nx::Int, ny::Int, nit::Int, rho::Float64)
        for n = 1:nit
                local pn = copy(p)
                for j = 2:ny-1
                        for i = 2:nx-1
                                p[i, j] = (
                                        (((pn[i+1, j] + pn[i-1, j]) * dy^2 + (pn[i, j+1] + pn[i, j-1]) * dx^2) / (2 * (dx^2 + dy^2)))
                                        - dx^2 * dy^2 / (2 * (dx^2 + dy^2))
                                        * rho * (1/dt * (
                                                (u[i+1, j] - u[i-1, j]) / (2 * dx) + (v[i, j+1] - v[i, j-1]) / (2 * dy)
                                                )
                                                - (((u[i+1, j] - u[i-1, j]) / (2 * dx)) ^ 2)
                                                - 2 * (u[i, j+1] - u[i, j-1]) / (2 * dy) * (v[i+1, j] - v[i-1, j]) / (2 * dx)
                                                - (((v[i, j+1] - v[i, j-1]) / (2 * dy)) ^ 2)
                                                )
                                        )
                        end
                end

                p[end, :] = p[end - 1, :] # dp/dx = 0 at x = 2
                p[:, 1] = p[:, 2] # dp/dy = 0 at y = 0
                p[1, :] = p[2, :] # dp/dx = 0 at x = 0
                p[:, end] .= 0 # p = 0 at y = 2
        end
        return nothing
end

# Note: I should make a mutable structure to contain all the config options. Would make my life a lot easier than having to input them all into every function.
function cavityFlow(nt::Int, u, v, p, rho::Float64, nu::Float64, nx::Int, ny::Int, dx::Float64, dy::Float64, nit::Int)
        for t = 1:nt
                removeDivergenceFromPressure!(p, u, v, dx, dy, dt, nx, ny, nit, rho)
                moveFluid!(u, v, p, nx, ny, dx, dy, dt, rho, nu)
        end
        return nothing
end


# Setting up some constant
nx = 41 # Number of steps in x
ny = 41 # Number of steps in y
nt = 300 # Number of steps in t
nit = 50 # Number of poisson iterations.

dx = 2 / (nx - 1) # Spatial step size
dy = 2 / (ny - 1)

x = range(0, 2, length = nx)
y = range(0, 2, length = ny)

rho = 1.0 # Density
nu = 0.1 # Viscosity
dt = 0.001 # Time step size

u = zeros(nx, ny)
v = zeros(nx, ny)
p = zeros(nx, ny)

cavityFlow(nt, u, v, p, rho, nu, nx, ny, dx, dy, nit)

heatmap(x, y, p)
# contour!(x, y, p)
#
arrows!(x, y, u, v, arrowsize = 0.02, lengthscale = 0.3)
