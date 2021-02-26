# Basic equations
# Velocity field
# Dv/Dt = -1/ρ ∇p + ν∇^2v + F (i.e. Navier-Stokes)

using LinearAlgebra
# using AbstractPlotting
# using GLMakie
using Plots

function moveFluid!(u, p, Fx, Fy, nx::Int, ny::Int, dx::Float64, dy::Float64, dt::Float64, rho, nu::Float64, periodicX::Bool, periodicY::Bool)
    local un = copy(u) # Make temporary copy of current velocity.
    local vn = copy(v) # Make these copies before calling the function.

    # Discretised advection / diffusion
    for j = 2:ny - 1
        for i = 2:nx - 1
            u[j, i] = (un[j, i]
                    - un[j, i] * (dt / dx) * (un[j, i] - un[j, i - 1])
                    - vn[j, i] * (dt / dy) * (un[j, i] - un[j - 1, i])
                    - dt / (rho[j, i] * 2 * dx) * (p[j, i+1] - p[j, i-1])
                    + nu * (
                        dt / (dx^2) * (un[j, i+1] - 2 * un[j, i] + un[j, i-1])
                        + dt / (dy^2) * (un[j+1, i] - 2 * un[j, i] + un[j-1, i])
                    )
                    + dt * Fx
                    )
            v[j, i] = (vn[j, i]
                    - un[j, i] * (dt / dx) * (vn[j, i] - vn[j, i - 1])
                    - vn[j, i] * (dt / dy) * (vn[j, i] - vn[j - 1, i])
                    - dt / (rho[j, i] * 2 * dy) * (p[j+1, i] - p[j-1, i])
                    + nu * (
                        dt / (dx^2) * (vn[j, i+1] - 2 * vn[j, i] + vn[j, i-1])
                        + dt / (dy^2) * (vn[j+1, i] - 2 * vn[j, i] + vn[j-1, i])
                    )
                    + dt * Fy
                    )
        end
    end

    velocityPeriodicBoundary!(u, un, vn, p, periodicX, periodicY)

    if !periodicX && !periodicY
        u[:, 1] .= -u[:, 2]
        u[:, end] .= -u[:, end-1]
        u[1, :] .= -u[2, :]
        u[end, :] .= -u[end-1, :]
        v[:, 1] .= -v[:, 2]
        v[:, end] .= -v[:, end-1]
        v[1, :] .= -v[2, :]
        v[end, :] .= -v[end-1, :]
    elseif periodicX && !periodicY
        u[1, :] .= 0
        u[end, :] .= 0
        v[1, :] .= 0
        v[end, :] .= 0
    elseif !periodicX && periodicY
        u[:, 1] .= 0
        u[:, end] .= 0
        v[:, 1] .= 0
        v[:, end] .= 0
    end

    return nothing
end

function removeDivergenceFromPressure!(p, u, v, nx::Int, ny::Int, dx::Float64, dy::Float64, dt::Float64, rho, targetError::Float64, periodicX::Bool, periodicY::Bool)
    currentError = 1
    # Iterate while currentError > targetError to converge pressure and remove divergence
    while currentError > targetError
    # for lk = 1:100
        local pn = copy(p)
        for j = 2:ny - 1
            for i = 2:nx - 1
                p[j, i] = (
                    ((pn[j, i+1] + pn[j, i-1]) * dy^2 + (pn[j+1, i] + pn[j-1, i]) * dx^2) / (2 * (dx^2 + dy^2))
                    - dx^2 * dy^2 / (2 * (dx^2 + dy^2))
                    * (rho[j, i]) * (1/dt * (
                        (u[j, i+1] - u[j, i-1]) / (2 * dx) + (v[j+1, i] - v[j-1, i]) / (2 * dy)
                        )
                        - ((u[j, i+1] - u[j, i-1]) / ((2 * dx)))^2
                        - (2 * (u[j+1, i] - u[j-1, i]) / (2 * dy)) * ((v[j, i+1] - v[j, i-1]) / (2*dx))
                        - ((v[j+1, i] - v[j-1, i]) / (2 * dy))^2
                    )
                )
            end
        end
        pressurePeriodicBoundary!(p, pn, u, v, periodicX, periodicY)
        # Remember to include non-periodic boundaries
        if !periodicX && !periodicY
            # Cavity boundaries
            p[:, end] .= p[:, end-1] - rho * nu / dx * (-2 * u[:, end-1] + u[:, end-2])
            p[:, 1] .= p[:, 2] - rho * nu / dx * (-2 * u[:, 2] + u[:, 3])
            p[1, :] .= p[2, :] - rho * nu / dy * (-2 * v[2, :] + v[3, :])
            p[end, :] .= p[end-1, :] - p[end-1, :] - rho * nu / dy * (-2 * u[end-1, :] + u[end-2, :])
        elseif periodicX && !periodicY
            p[1, :] .= p[2, :] - rho * nu / dy * (-2 * v[2, :] + v[3, :])
            p[end, :] .= p[end-1, :] - rho * nu / dy * (-2 * v[end-1, :] + v[end-2, :])
        elseif !periodicX && periodicY
            p[:, end] .= p[:, end-1] - rho * nu / dx * (-2 * u[:, end-1] + u[:, end-2])
            p[:, 1] .= p[:, 2] - rho * nu / dx * (-2 * u[:, 2] + u[:, 3])
        end

        currentError = (sum(abs.(p)) - sum(abs.(pn))) / sum(abs.(p))
    end
    return nothing
end

function velocityPeriodicBoundary!(u, un, vn, p, horizontal::Bool, vertical::Bool)
    # For horizontal boundary conditions, solve for i=1, i=nx, j=2:ny-1
    # For vertical boundary conditions, solve for j=1, j=ny, i=2:nx-1
    # How do to corners????

    if horizontal
        for j = 2:ny - 1
            for i = (1, nx)
                u[j, i] = (un[j, i]
                        - un[j, i] * (dt / dx) * (un[j, i] - un[j, mod1(i-1, nx)])
                        - vn[j, i] * (dt / dy) * (un[j, i] - un[j - 1, i])
                        - dt / (rho[j, i] * 2 * dx) * (p[j, mod1(i+1, nx)] - p[j, mod1(i-1, nx)])
                        + nu * (
                            dt / (dx^2) * (un[j, mod1(i+1, nx)] - 2 * un[j, i] + un[j, mod1(i-1, nx)])
                            + dt / (dy^2) * (un[j+1, i] - 2 * un[j, i] + un[j-1, i])
                        )
                        + dt * Fx
                        )
                v[j, i] = (vn[j, i]
                        - un[j, i] * (dt / dx) * (vn[j, i] - vn[j, mod1(i-1, nx)])
                        - vn[j, i] * (dt / dy) * (vn[j, i] - vn[j - 1, i])
                        - dt / (rho[j, i] * 2 * dy) * (p[j+1, i] - p[j-1, i])
                        + nu * (
                            dt / (dx^2) * (vn[j, mod1(i+1, nx)] - 2 * vn[j, i] + vn[j, mod1(i-1, nx)])
                            + dt / (dy^2) * (vn[j+1, i] - 2 * vn[j, i] + vn[j-1, i])
                        )
                        + dt * Fy
                        )
            end
        end
    end
    if vertical
        for j = (1, ny)
            for i = 2:nx - 1
                u[j, i] = (un[j, i]
                        - un[j, i] * (dt / dx) * (un[j, i] - un[j, i-1])
                        - vn[j, i] * (dt / dy) * (un[j, i] - un[mod1(j-1, ny), i])
                        - dt / (rho[j, i] * 2 * dx) * (p[j, i+1] - p[j, i-1])
                        + nu * (
                            dt / (dx^2) * (un[j, i+1] - 2 * un[j, i] + un[j, i-1])
                            + dt / (dy^2) * (un[mod1(j+1, ny), i] - 2 * un[j, i] + un[mod1(j-1, ny), i])
                        )
                        + dt * Fx
                        )
                v[j, i] = (vn[j, i]
                        - un[j, i] * (dt / dx) * (vn[j, i] - vn[j, i-1])
                        - vn[j, i] * (dt / dy) * (vn[j, i] - vn[mod1(j-1, ny), i])
                        - dt / (rho[j, i] * 2 * dy) * (p[mod1(j+1, ny), i] - p[mod1(j-1, ny), i])
                        + nu * (
                            dt / (dx^2) * (vn[j, i+1] - 2 * vn[j, i] + vn[j, i-1])
                            + dt / (dy^2) * (vn[mod1(j+1, ny), i] - 2 * vn[j, i] + vn[mod1(j-1, ny), i])
                        )
                        + dt * Fy
                        )
            end
        end
    end
    if horizontal && vertical
        # set corners.
        for j = (1, ny)
            for i = (1, nx)
                u[j, i] = (un[j, i]
                        - un[j, i] * (dt / dx) * (un[j, i] - un[j, mod1(i-1, nx)])
                        - vn[j, i] * (dt / dy) * (un[j, i] - un[mod1(j-1, ny), i])
                        - dt / (rho[j, i] * 2 * dx) * (p[j, mod1(i+1, nx)] - p[j, mod1(i-1, nx)])
                        + nu * (
                            dt / (dx^2) * (un[j, mod1(i+1, nx)] - 2 * un[j, i] + un[j, mod1(i-1, nx)])
                            + dt / (dy^2) * (un[mod1(j+1, ny), i] - 2 * un[j, i] + un[mod1(j-1, ny), i])
                        )
                        + dt * Fx
                        )
                v[j, i] = (vn[j, i]
                        - un[j, i] * (dt / dx) * (vn[j, i] - vn[j, mod1(i-1, nx)])
                        - vn[j, i] * (dt / dy) * (vn[j, i] - vn[mod1(j-1, ny), i])
                        - dt / (rho[j, i] * 2 * dy) * (p[mod1(j+1, ny), i] - p[mod1(j-1, ny), i])
                        + nu * (
                            dt / (dx^2) * (vn[j, mod1(i+1, nx)] - 2 * vn[j, i] + vn[j, mod1(i-1, nx)])
                            + dt / (dy^2) * (vn[mod1(j+1, ny), i] - 2 * vn[j, i] + vn[mod1(j-1, ny), i])
                        )
                        + dt * Fy
                        )
            end
        end
    end

end

function pressurePeriodicBoundary!(p, pn, u, v, horizontal::Bool, vertical::Bool)
    if horizontal
        for j = 2:ny - 1
            for i = (1, nx)
                p[j, i] = (
                    ((pn[j, mod1(i+1, nx)] + pn[j, mod1(i-1, nx)]) * dy^2 + (pn[j+1, i] + pn[j-1, i]) * dx^2) / (2 * (dx^2 + dy^2))
                    - dx^2 * dy^2 / (2 * (dx^2 + dy^2))
                    * (rho[j, i]) * (1/dt * (
                        (u[j, mod1(i+1, nx)] - u[j, mod1(i-1, nx)]) / (2 * dx) + (v[j+1, i] - v[j-1, i]) / (2 * dy)
                        )
                        - ((u[j, mod1(i+1, nx)] - u[j, mod1(i-1, nx)]) / ((2 * dx)))^2
                        - (2 * (u[j+1, i] - u[j-1, i]) / (2 * dy)) * ((v[j, mod1(i+1, nx)] - v[j, mod1(i-1, nx)]) / (2*dx))
                        - ((v[j+1, i] - v[j-1, i]) / (2 * dy))^2
                    )
                )
            end
        end
    end
    if vertical
        for j = (1, ny)
            for i = 2:nx - 1
                p[j, i] = (
                    ((pn[j, i+1] + pn[j, i-1]) * dy^2 + (pn[mod1(j+1, ny), i] + pn[mod1(j-1, ny), i]) * dx^2) / (2 * (dx^2 + dy^2))
                    - dx^2 * dy^2 / (2 * (dx^2 + dy^2))
                    * (rho[j, i]) * (1/dt * (
                        (u[j, i+1] - u[j, i-1]) / (2 * dx) + (v[mod1(j+1, ny), i] - v[mod1(j-1, ny), i]) / (2 * dy)
                        )
                        - ((u[j, i+1] - u[j, i-1]) / ((2 * dx)))^2
                        - (2 * (u[mod1(j+1, ny), i] - u[mod1(j-1, ny), i]) / (2 * dy)) * ((v[j, i+1] - v[j, i-1]) / (2*dx))
                        - ((v[mod1(j+1, ny), i] - v[mod1(j-1, ny), i]) / (2 * dy))^2
                    )
                )
            end
        end
    end
    if horizontal && vertical
        # set corners.
        for j = (1, ny)
            for i = (1, nx)
                p[j, i] = (
                    ((pn[j, mod1(i+1, nx)] + pn[j, mod1(i-1, nx)]) * dy^2 + (pn[mod1(j+1, ny), i] + pn[mod1(j-1, ny), i]) * dx^2) / (2 * (dx^2 + dy^2))
                    - dx^2 * dy^2 / (2 * (dx^2 + dy^2))
                    * (rho[j, i]) * (1/dt * (
                        (u[j, mod1(i+1, nx)] - u[j, mod1(i-1, nx)]) / (2 * dx) + (v[mod1(j+1, ny), i] - v[mod1(j-1, ny), i]) / (2 * dy)
                        )
                        - ((u[j, mod1(i+1, nx)] - u[j, mod1(i-1, nx)]) / ((2 * dx)))^2
                        - (2 * (u[mod1(j+1, ny), i] - u[mod1(j-1, ny), i]) / (2 * dy)) * ((v[j, mod1(i+1, nx)] - v[j, mod1(i-1, nx)]) / (2*dx))
                        - ((v[mod1(j+1, ny), i] - v[mod1(j-1, ny), i]) / (2 * dy))^2
                    )
                )
            end
        end
    end
end

function moveDensity!(rho, u, v, dx, dy, dt, nx, ny)
    local rhon = copy(rho)

    for j = 2:ny - 1
        for i = 2:nx - 1
            rho[j, i] = (rhon[j, i]
                        - dt/(2*dx) * (u[j, i] * (rhon[j, i+1] + rhon[j, i]) - u[j, i-1] * (rhon[j, i] + rhon[j, i-1]))
                        - dt/(2*dy) * (v[j, i] * (rhon[j+1, i] + rhon[j, i]) - u[j-1, i] * (rhon[j, i] + rhon[j-1, i]))
                        + dt/(dx^2) * (rhon[j, i+1] - 2 * rhon[j, i] + rhon[j, i-1])
                        + dt/(dy^2) * (rhon[j+1, i] - 2 * rhon[j, i] + rhon[j-1, i])
                        )

            # (rhon[j, i]
            #             - dt/dx * u[j, i] * (rhon[j, i] - rhon[j, i-1])
            #             - dt/dy * v[j, i] * (rhon[j, i] - rhon[j, i-1])
            #             - dt/dx * rhon[j, i] * (u[j, i] - u[j, i-1])
            #             - dt/dy * rhon[j, i] * (v[j, i] - v[j-1, i])
            #             )
        end
    end
    return nothing
end

function moveDatFluidBro(p, u, v, nx, ny, dx, dy, dt, rho, targetError, Fx, Fy, nu, nt, periodicX::Bool, periodicY::Bool)
    for t = 1:nt
        # Update pressure
        removeDivergenceFromPressure!(p, u, v, nx, ny, dx, dy, dt, rho, targetError, periodicX, periodicY)
        # Update velocities
        moveFluid!(u, p, Fx, Fy, nx, ny, dx, dy, dt, rho, nu, periodicX, periodicY)
        # Update density
        moveDensity!(rho, u, v, dx, dy, dt, nx, ny)
        percent = t/nt*100
        println(string(percent, '%'))
    end
    return nothing
end

## Set up initial states, grids, define terms, etc

nx = 15 # Number of x grid points
ny = 15 # Number of y grid points
nt = 10 # Number of time points
targetError = 1e-8 # Target error for the pressure divergence removal

dx = 2 / (nx - 1) # Step size in x
dy = 2 / (ny - 1) # Step size in y

x = range(0, 1, length = nx)
y = range(0, 1, length = ny)

rho1 = 1.0 # Density of not drop
rho2 = 2.0 # Density of drop
nu = 0.02 # Viscosity
dt = 0.00125 # Time step size
Fx = 0
Fy = -100

# Set periodic boundaries on x/y
periodicX = true
periodicY = true
# Both false defaults to a simple cavity
# Otherwise, defaults to channel in X/Y direction

u = zeros(ny, nx) # Horizontal velocity
v = zeros(ny, nx) # Vertical velocity
p = zeros(ny, nx) # Pressure
rho = zeros(ny, nx) .+ rho1 # Density of everything

for j = 2:ny-1
    for i = 2:nx-1
        if (i - 7)^2 + (j - 7)^2 <= 2^2
            rho[j, i] = rho2 # Set density of drop
        end
    end
end

moveDatFluidBro(p, u, v, nx, ny, dx, dy, dt, rho, targetError, Fx, Fy, nu, nt, periodicX, periodicY)

heatmap(x, y, p)
