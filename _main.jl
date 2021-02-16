using LinearAlgebra
using AbstractPlotting
using GLMakie
# AbstractPlotting.inline!(false)

# Basic equations
# Velocity field
# Dv/Dt = -1/ρ ∇p + ν∇^2v + F (i.e. Navier-Stokes)
# Density
# Dρ/Dt + ρ∇⋅v = 0
# Composition
# DS/Dt = S{dotted}

# Onto NS with a "source" term. i.e. a force field.

function applyPeriodicVelXBoundary!(a, an, bn, p, F, i, nx::Int, ny::Int)
        for j = 2:ny-1
                a[j, i] = (an[j, i] # Discretised Du/Dt = -1/ρ ∇p + ν∇^2u
                        - an[j, i] * (dt / dx) * (an[j, i] - an[j, mod1(i-1, nx)])
                        - bn[j, i] * (dt / dy) * (an[j, i] - an[j-1, j])
                        - dt / (rho * 2 * dx) * (p[j, mod1(i+1, nx)] - p[j, mod1(i-1, nx)])
                        + nu * (
                                dt / (dx^2) * (an[j, mod1(i+1, nx)] - 2 * an[j, i] + an[j, mod1(i-1, nx)])
                                + dt / (dy^2) * (an[j+1, i] - 2 * an[j, i] + an[j-1, i])
                                )
                        + dt * F # Force only in horizontal direction, since we're looking at a horizontal channel.
                        )
        end
        return nothing
end

function applyPeriodicPXBoundary!(p, pn, u, v, i, nx::Int, ny::Int)
        for j = 2:ny-1
                p[j, i] = (
                        (((pn[j, mod1(i+1, nx)] + pn[j, mod1(i-1, nx)]) * dy^2 + (pn[j+1, i] + pn[j-1, i]) * dx^2) / (2 * (dx^2 + dy^2)))
                        - dx^2 * dy^2 / (2 * (dx^2 + dy^2))
                        * rho * (1/dt * (
                                (u[j, mod1(i+1, nx)] - u[j, mod1(i-1, nx)]) / (2 * dx) + (v[j+1, i] - v[j-1, i]) / (2 * dy)
                                )
                                - (((u[j, mod1(i+1, nx)] - u[j, mod1(i-1, nx)]) / (2 * dx)) ^ 2)
                                - 2 * (u[j+1, i] - u[j-1, i]) / (2 * dy) * (v[j, mod1(i+1, nx)] - v[j, mod1(i-1, nx)]) / (2 * dx)
                                - (((v[j+1, i] - v[j-1, i]) / (2 * dy)) ^ 2)
                                )
                        )
        end
        return nothing # Libbum tells me that this is faster than actually not returning anything.
end

function moveFluid!(u, v, p, F, nx::Int, ny::Int, dx::Float64, dy::Float64, dt::Float64, rho::Float64, nu::Float64)
        local un = copy(u) # Make temporary copy of previous time step
        local vn = copy(v)
        for j = 2:ny - 1
                for i = 2:nx - 1
                        u[j, i] = (un[j, i] # Discretised Du/Dt = -1/ρ ∇p + ν∇^2u
                                - un[j, i] * (dt / dx) * (un[j, i] - un[j, i-1])
                                - vn[j, i] * (dt / dy) * (un[j, i] - un[j-1, j])
                                - dt / (rho * 2 * dx) * (p[j, i+1] - p[j, i-1])
                                + nu * (
                                        dt / (dx^2) * (un[j, i+1] - 2 * un[j, i] + un[j, i-1])
                                        + dt / (dy^2) * (un[j+1, i] - 2 * un[j, i] + un[j-1, i])
                                        )
                                + dt * F # Force only in horizontal direction, since we're looking at a horizontal channel.
                                )

                        v[j, i] = (vn[j, i] # Discretised Dv/Dt = -1/ρ ∇p + ν∇^2v
                                - un[j, i] * (dt / dx) * (vn[j, i] - vn[j, i-1])
                                - vn[j, i] * (dt / dy) * (vn[j, i] - vn[j-1, j])
                                - dt / (rho * 2 * dy) * (p[j+1, i] - p[j-1, i])
                                + nu * (
                                        dt / (dx^2) * (vn[j, i+1] - 2 * vn[j, i] + vn[j, i-1])
                                        + dt / (dy^2) * (vn[j+1, i] - 2 * vn[j, i] + vn[j-1, i])
                                        )
                                )

                                # Probably actually can put the un / vn in as an input. Yeah, that probably would have been the smarter thing to do
                                # Especially when this gets to 3D.
                end
        end

        # Left and right side *are* periodic, and need to be treated carefully. A naiive setting of the values at one side to the other might work, but an explicit calculation is better.
        j = 1
        applyPeriodicVelXBoundary!(u, un, vn, p, F, j, nx, ny)
        applyPeriodicVelXBoundary!(v, vn, un, p, 0, j, nx, ny)
        j = ny
        applyPeriodicVelXBoundary!(u, un, vn, p, F, j, nx, ny)
        applyPeriodicVelXBoundary!(v, vn, un, p, 0, j, nx, ny)

        # Apply boundary conditions.
        # Bottom and top are not periodic, and are 0.
        u[1, :] .= 0
        u[end, :] .= 0
        v[1, :] .= 0
        v[end, :] .= 0

        return un
end

# Poisson should take in p, u, v, dx, dy, dt. The version in the tutorial has a set number of iterations (nit), so maybe the err doesn't converge nicely?
# This is also where I will eventually put Jesse's nice stencils. Also just realised that iterative to equilibrium removes the divergence! Of course it does!
# I can be really thick sometimes.

function removeDivergenceFromPressure!(p, u, v, dx::Float64, dy::Float64, dt::Float64, nx::Int, ny::Int, nit::Int, rho::Float64)
        for n = 1:nit
                local pn = copy(p)
                for j = 2:ny-1
                        for i = 2:nx-1
                                p[j, i] = (
                                        (((pn[j, i+1] + pn[j, i-1]) * dy^2 + (pn[j+1, i] + pn[j-1, i]) * dx^2) / (2 * (dx^2 + dy^2)))
                                        - dx^2 * dy^2 / (2 * (dx^2 + dy^2))
                                        * rho * (1/dt * (
                                                (u[j, i+1] - u[j, i-1]) / (2 * dx) + (v[j+1, i] - v[j-1, i]) / (2 * dy)
                                                )
                                                - (((u[j, i+1] - u[j, i-1]) / (2 * dx)) ^ 2)
                                                - 2 * (u[j+1, i] - u[j-1, i]) / (2 * dy) * (v[j, i+1] - v[j, i-1]) / (2 * dx)
                                                - (((v[j+1, i] - v[j-1, i]) / (2 * dy)) ^ 2)
                                                )
                                        )
                        end
                end
                # Periodic boundaries at left and right
                j = 1
                applyPeriodicPXBoundary!(p, pn, u, v, j, nx, ny)
                j = ny
                applyPeriodicPXBoundary!(p, pn, u, v, j, nx, ny)

                # Flux is 0 at top and bottom
                p[1, :] = p[2, :] # dp/dy = 0 at y = 0
                p[end, :] = p[end-1, :] # p = 0 at y = 2

        end
        return nothing
end

# Note: I should make a mutable structure to contain all the config options. Would make my life a lot easier than having to input them all into every function.
function cavityFlow(err::Float64, u, v, p, F, rho::Float64, nu::Float64, nx::Int, ny::Int, dx::Float64, dy::Float64, nit::Int)
        currentErr = 1
        while currentErr > err
                removeDivergenceFromPressure!(p, u, v, dx, dy, dt, nx, ny, nit, rho)
                un = moveFluid!(u, v, p, F, nx, ny, dx, dy, dt, rho, nu)
                currentErr = (sum(u) - sum(un)) / sum(u)
        end
        return numsteps
end


# Setting up some constant
nx = 41 # Number of steps in x
ny = 41 # Number of steps in y
err = 0.001 # Target error
nit = 50 # Number of poisson iterations.

dx = 2 / (nx - 1) # Spatial step size
dy = 2 / (ny - 1)

x = range(0, 2, length = nx)
y = range(0, 2, length = ny)

rho = 1.0 # Density
nu = 0.1 # Viscosity
dt = 0.001 # Time step size
F = 100

u = zeros(ny, nx)
v = zeros(ny, nx)
p = zeros(ny, nx)

cavityFlow(err, u, v, p, F, rho, nu, nx, ny, dx, dy, nit)

contour(y, x, u') # Transposed so it matches how we want it to plot - horizontally.
