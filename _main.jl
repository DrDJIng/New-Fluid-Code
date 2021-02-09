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
dx = 2 * pi / (nx - 1)
nt = 500 # Number of timesteps
nu = 0.07 # Viscosity
sigma = 0.2 # Sigma is randomly defined later, then we're told to ignore it for now. K. Probably some kind of discretised stepping rate.
dt = dx * nu # dt is the amount of time each step, i.e. delta t

t = 0 # Initial time
x = range(0, 2*pi, length = nx)

# Set up initial conditions. We want 1 everywhere, except where there are twos (0.5 ≤ x ≤ 1)
u = ones(1, nx)
u[round(Int, 0.5/dx):round(Int, 1/dx+1)] .= 2

phi = exp.(-x.*x./(4*nu)) + exp.(-(x .- 2 * pi).^2 ./ (4 * nu))
dphidx = -(-8*t .+ 2*x) .* exp.(-(-4*t .+ x).^2 ./ (4*nu*(t + 1))) ./ (4*nu*(t + 1)) .- (-8*t .+ 2*x .- 4*pi) .* exp.(-(-4*t .+ x .- 2*pi).^2 ./ (4*nu .* (t + 1))) ./ (4*nu*(t + 1))
u = -2*nu.*(-(-8*t .+ 2*x).*exp.(-(-4*t .+ x).^2 ./ (4*nu*(t + 1))) ./ (4*nu*(t + 1)) .- (-8*t .+ 2*x .- 4*pi) .* exp.(-(-4*t .+ x .- 2*pi).^2 ./ (4*nu .* (t + 1))) ./ (4*nu .* (t + 1))) ./ (exp.(-(-4*t .+ x .- 2*pi) .^2 ./ (4*nu * (t + 1))) .+ exp.(-(-4*t .+ x).^2 ./ (4*nu*(t + 1)))) .+ 4
# Copied directly from the tutorial and modified for Julia, because screeeeewwwww typing this out by hand.

uini = copy(u)
plot(uini)

# For every element in u, we need to perform the discretization
anim = @animate for n in 1:nt
# for n in 1:nt # This n actually isn't doing anything at the moment. I guess we *could* use it to store time information, but atm we just lose that.
    local un = copy(u)
    for i in 2:nx-1
        u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + nu * dt / (dx * dx) * (un[i+1] - 2 * un[i] + un[i-1]) # Discretised diffusion added here. It's linear, but based on viscosity.
        u[1] = un[1] - un[1] * dt / dx * (un[1] - un[end - 1]) + nu * dt / (dx * dx) * (un[2] - 2 * un[1] + un[end - 1])
        u[end] = u[1]
    end
    plot(uini)
    plot!(u)
end

# plot!(u)
gif(anim, "sawtooth.gif", fps = 30)
