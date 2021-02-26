# using Plots
using GLMakie

# Domain size and physical variables
Lx = 1.0
Ly = 1.0
gx = 0.0
gy = -100.0
rho1 = 1.0
rho2 = 2.0

# Dynamic viscosity
m0 = 0.05

# Tangential velocities
uNorth = 0.0
uSouth = 0.0
vWest = 0.0
vEast = 0.0

# Numerical variables
nx = 21
ny = 21
dt = 4.0e-4
nstep = 400
maxiter = 1000
maxError = 1.0e-3
beta = 1.2

u = zeros(nx+2, ny+3)
v = zeros(nx+3, ny+2)
p = zeros(nx+3, ny+3)

ut = zeros(nx+2, ny+3)
vt = zeros(nx+3, ny+2)

tmp1 = zeros(nx+2, ny+2)
tmp2 = zeros(nx+2, ny+2)

# Velocities at center of grid for plotting
uu = zeros(nx+2, ny+2)
vv = zeros(nx+2, ny+2)
x = range(0, step = Lx, stop = nx)
y = range(0, step = Ly, stop = ny)

# Define the grid
dx = Lx / nx
dy = Ly / ny

r = zeros(nx+3, ny+3) .+ rho1

for j = 1:ny+3
    for i = 1:nx+3
        if (i - 12)^2 + (j - 12)^2 <= 2^2
            r[i, j] = rho2 # Set density of drop
        end
    end
end

time = 0.0

for steps in range(1, stop = nstep)

    # Set the tangential velocities at the boundaries.
    u[:, 1] = 2 * uSouth .- u[:, 2]
    u[:, end] = 2 * uNorth .- u[:, end-1]
    v[1, :] = 2 * vWest .- v[1, :]
    v[end, :] = 2 * vEast .- v[end-1, :]

    # Temporary velocities
    # U-velocity
    for i in range(2, stop = nx+1)
        for j in range(2, stop = ny+2)
            ut[i,j] = u[i,j]+dt*(-0.25*(((u[i+1,j]+u[i,j])^2-(u[i,j]+
                        u[i-1,j])^2)/dx+((u[i,j+1]+u[i,j])*(v[i+1,j]+
                        v[i,j])-(u[i,j]+u[i,j-1])*(v[i+1,j-1]+v[i,j-1]))/dy)+
                        m0/(0.5*(r[i+1,j]+r[i,j]))*((u[i+1,j]-2*u[i,j]+u[i-1,j])/dx^2+
                        (u[i,j+1]-2*u[i,j]+u[i,j-1])/dy^2 )+gx)
        end
    end

    for i in range(2, stop = nx+2)
        for j in range(2, stop = ny+1)
            vt[i,j] = v[i,j]+dt*(-0.25*(((u[i,j+1]+u[i,j])*(v[i+1,j]+
                v[i,j])-(u[i-1,j+1]+u[i-1,j])*(v[i,j]+v[i-1,j]))/dx+
                ((v[i,j+1]+v[i,j])^2-(v[i,j]+v[i,j-1])^2)/dy)+
                m0/(0.5*(r[i,j+1]+r[i,j]))*((v[i+1,j]-2*v[i,j]+v[i-1,j])/dx^2+
                (v[i,j+1]-2*v[i,j]+v[i,j-1])/dy^2 )+gy)
        end
    end

    # Compute source term and coefficient for p[i, j]
    rt = copy(r)
    lrg = 1000.0
    rt[:, 1] .= lrg
    rt[:, end] .= lrg
    rt[1, :] .= lrg
    rt[end, :] .= lrg

    for i in range(2, stop = nx+2)
        for j in range(2, stop = ny+2)
            tmp1[i,j] = (0.5/dt)*( (ut[i,j]-ut[i-1,j])/dx+(vt[i,j]-vt[i,j-1])/dy)
            tmp2[i,j] = 1.0/( (1/dx)*( 1/(dx*(rt[i+1,j]+rt[i,j]))+
                1/(dx*(rt[i-1,j]+rt[i,j])))+
                (1/dy)*(1/(dy*(rt[i,j+1]+rt[i,j]))+
                1/(dy*(rt[i,j-1]+rt[i,j]))))
        end
    end

    iter = 0
    while true
        pn = copy(p)
        iter = iter+1
        for i in range(2, stop = nx+2)
            for j in range(2, stop = ny+2)
                p[i,j]=(1.0-beta)*p[i,j]+beta*tmp2[i,j]*(
                    (1/dx)*( p[i+1,j]/(dx*(rt[i+1,j]+rt[i,j]))+
                    p[i-1,j]/(dx*(rt[i-1,j]+rt[i,j])))+
                    (1/dy)*( p[i,j+1]/(dy*(rt[i,j+1]+rt[i,j]))+
                    p[i,j-1]/(dy*(rt[i,j-1]+rt[i,j])))-tmp1[i,j])
            end
        end
        if maximum(abs.(pn.-p))<maxError
            break
        end
        if iter>maxiter
            break
        end
    end

    #CORRECT THE u-velocity
    for i in range(2, stop = nx)
        for j in range(2, stop = ny+2)
            u[i,j] = ut[i,j]-dt*(2.0/dx)*(p[i+1,j]-p[i,j])/(r[i+1,j]+r[i,j])
        end
    end

    #CORRECT THE v-velocity
    for i in range(2, stop = nx+2)
        for j in range(2, stop = ny)
            v[i,j] = vt[i,j]-dt*(2.0/dy)*(p[i,j+1]-p[i,j])/(r[i,j+1]+r[i,j])
        end
    end

# ADVECT DENSITY using centered difference plus diffusion
    ro = copy(r)
    for i in range(2, stop = nx+2)
        for j in range(2, stop = ny+2)
            r[i,j] = ro[i,j]-(0.5*dt/dx)*(u[i,j]*(ro[i+1,j]
                    +ro[i,j])-u[i-1,j]*(ro[i-1,j]+ro[i,j]))
                    -(0.5* dt/dy)*(v[i,j]*(ro[i,j+1]
                    +ro[i,j])-v[i,j-1]*(ro[i,j-1]+ro[i,j]))
                    +(m0*dt/dx/dx)*(ro[i+1,j]-2.0*ro[i,j]+ro[i-1,j])
                    +(m0*dt/dy/dy)*(ro[i,j+1]-2.0*ro[i,j]+ro[i,j-1])
        end
    end
    time = time+dt

end



uu = 0.5 * (u[1:nx+1, 2:ny+2] + u[1:nx+1, 1:ny+1])
vv = 0.5 * (v[2:nx+2, 1:ny+1] + v[1:nx+1, 1:ny+1])
fig = contour(x, y, r[2:end-1, 2:end-1]) #log.(sqrt.(uu.^2 + vv.^2))
arrows!(x, y, uu, vv, lengthscale = 1.0)
fig
