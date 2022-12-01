###################################################################
#                                                                 #
#         Computational Fluid Dynamics (24-718) - Methods         #
#                                                                 #
###################################################################


import numpy as np


def SolveAdvectionDiffusion(T, u, v, dx, dy, Nx, Ny, q, k, rho, Cp):
    # Pouya
    T_t = T.copy()

    #Interior Nodes
    T_t[2:Nx-1,2:Ny-1] = 

    #Boundary Nodes (with Neumann BC everywhere except at top)
    #Left and Right Boundary
    T_t[1,1] =
    T_t[Nx,1] =
    for j in range(2:Ny-1):
        T_t[1,j] =
        T_t[Nx-1,j] =
    T_t[1,Ny-1] =
    T_t[Nx,Ny-1] =

    # Top and Bottom
    for l in range(2:Nx-1):
        T_t[l,1] =
        T_t[l,Ny-1] =

    return T_t

def calPoisson(u_div, dt, rho, Nx, Ny, tol):
    p = np.zeros([Nx+1,Ny+1])
    err = 1
    # Jinho
    while err > tol:
        newp[1:Nx,1:Ny] = 
        err = np.max(np.abs(newp - p))
        p = newp.copy()
    return p

def calPressureGradient(p, dx, dy, Nx, Ny):
    # Pouya
    return dp_x, dp_y

def calVelocityDivergence(us, vs, dx, dy, Nx, Ny):
    # Pouya
    return u_div

def NavierStokesTemperature(u_bc_h, u_bc_v, v_bc_h, v_bc_v, Ta, nu, dt, dx, dy, lx, ly, \
            Nx, Ny, k, rho, Cp, num_timesteps, q, tol):

    # Define velocity field
    u = np.zeros([Nx+1,Ny+1])
    v = np.zeros([Nx+1,Ny+1])
    T = Ta*np.ones([Nx+1,Ny+1])
    
    # Initialize Dirichlet BCs
    u[:,0] = u_bc_v[:,0]
    u[:,-1] = u_bc_v[:,1]
    u[0,1:Ny] = u_bc_h[1:Ny,0]
    u[-1,1:Ny] = u_bc_h[1:Ny,1]

    for n in range(1,num_timesteps):
        # Calculate u* and v*
        us = u.copy()
        vs = v.copy()
        
        # For interior nodes (2:Nx-2), (2:Ny-2)
        for l in range(2,Nx-1):
            for j in range(2:Ny-1):
                us[l,j] = u[l,j] \
                        + (dt/(4*dx))*(u[l+1,j] - u[l-1,j])*(u[l+1,j] + u[l-1,j] + 2*u[l,j]) \
                        + (dt/(4*dy))*(u[l,j]*v[l,j] + u[l,j]*v[l+1,j] + u[l,j+1]*v[l,j] - u[l,j]*v[l,j-1] \
                                      -u[l,j]*v[l+1,j-1] - u[l,j-1]*v[l,j-1] - u[l,j-1]*v[l+1,j-1]) \
                        + 0.5*nu*(dt/dx**2)*(u[l+1,j] - 2*u[l,j] + u[l-1,j])  \
                        + 0.5*nu*(dt/dy**2)*(u[l,j+1] - 2*u[l,j] + u[l,j-1])
                vs[l,j] = 

        # Boundary adjacent nodes (1,Nx-1), (1,Ny-1)
        # Deal with Neumann and Dirichlet BCs here
        #Left and Right Boundary
        us[1,1] = 
        us[Nx,1] = 
        vs[1,1] = 
        vs[Nx,1] =
        for j in range(2:Ny-1):
            us[1,j] = 
            us[Nx-1,j] =
            vs[1,j] = 
            vs[Nx-1,j] = 
        us[1,Ny-1] =
        us[Nx,Ny-1] = 
        vs[1,Ny-1] = 
        vs[Nx,Ny-1] = 

        # Top and Bottom
        for l in range(2:Nx-1):
            us[l,1] = 
            us[l,Ny-1] = 
            vs[l,1] = 
            vs[l,Ny-1] =

        # Calculate Divergence of us,vs
        u_div = calVelocityDivergence(us, vs, dx, dy, Nx, Ny)

        # Calculate pressure by solving Poisson equation
        p = calPoisson(u_div, dt, rho, Nx, Ny, tol)

        # Calculate pressure gradient
        dp_x, dp_y = calPressureGradient(p, dx, dy, Nx, Ny)

        # Calculate new velocity
        u[1:Nx,1:Ny] = us[1:Nx,1:Ny] - (dt/rho)*dp_x[1:Nx,1:Ny]
        v[1:Nx,1:Ny] = vs[1:Nx,1:Ny] - (dt/rho)*dp_y[1:Nx,1:Ny]

        # Set Neumann BC values for velocity
        u[0,:int(d/dy)+1] = u[1,:int(d/dy)+1]
        u[Nx,:int(d/dy)+1] = u[Nx-1,:int(d/dy)+1]
        v[0,:int(d/dy)+1] = v[1,:int(d/dy)+1]
        v[Nx,:int(d/dy)+1] = v[Nx-1,:int(d/dy)+1]

        # Solve Advection-Diffusion equation for temperature
        T = SolveAdvectionDiffusion(T, u, v, dx, dy, Nx, Ny, q, k, rho, Cp)
    return u, v, T
