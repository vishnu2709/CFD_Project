#######################################################################
#                                                                     #
#          Computational Fluid Dynamics (24-718) - Project            #
#                                                                     #
#######################################################################


# Import required python modules
import numpy as np
import time

# Import additional modules
import methods
import plots


if __name__ == "__main__":
    # Parameters
    lx = 0.04 #m
    ly = 0.04 #m
    d = 0.01 #m
    dprime = 0.01 #m
    nu = 0.01 #m^2/s
    V = 2 # m/s
    Cp = 1005 #J/kg-K
    rho = 1.204 #kg/m^3
    k = 0.0262 # W/m-K
    q = 100000 # W/m^2
    t_max = 5 # s

    dx = 0.001 # m
    dy = 0.001 # m
    dt = 0.001 # s

    tol = 1e-5
    Nx = int(lx/dx)
    Ny = int(ly/dy)
    num_timesteps = int(t_max/dt) + 1

    # Initial Conditions
    uo = 0 # m/s
    vo = 0 # m/s
    To = 300 # K

    # Define boundary conditions
    # Velocity (Neumann parts will be dealt with later)

    u_bc_v = np.zeros([Nx+1,2])
    u_bc_h = np.zeros([Ny+1,2])

    v_bc_v = np.zeros([Nx+1,2])
    v_bc_h = np.zeros([Ny+1,2])


    # left boundary
    u_bc_h[:,0] = 0
    # right boundary
    u_bc_h[:,1] = 0
    # bottom boundary
    u_bc_v[:,0] = 0
    # top boundary
    u_bc_v[:,1] = 0

    # left boundary
    v_bc_h[:,0] = 0
    # right boundary
    v_bc_h[:,1] = 0
    # bottom boundary
    v_bc_v[:,0] = 0
    # top boundary
    v_bc_v[:,1] = 0
    v_bc_v[int((Lx-d)/(2*dx))+1:int((Lx+d)/(2*dx))+1,1] = V

    # Temperature (Dirichlet only)
    Ta = 300 #K


    ###############################################################

    
    print('Number of points (x-direction): {0:2d} '.format(Nx+1))
    print('Number of points (y-direction): {0:2d} '.format(Ny+1))
    print('Mesh size (dx): {0:.8f} mm'.format(dx))
    print('Mesh size (dy): {0:.8f} mm'.format(dy))
    print('Number of time steps: {0:2d} '.format(num_timesteps))
    print('Time step (dt): {0:.8f} s'.format(dt))


    start = time.time()
    u_t, v_t, T_t = methods.NavierStokesTemperature(u_bc_h, u_bc_v, v_bc_h, v_bc_v, Ta, nu, dt, dx, dy, lx, ly, \
            Nx, Ny, k, rho, Cp, num_timesteps, q, tol)
    end = time.time()
    print('time elpased: {0:.8f} s'.format(end - start))
    
    plots.VelocityFiled(np.linspace(0,lx,Nx+1),np.linspace(0,ly,Ny+1),u_t.transpose(),v_t.transpose())
    plots.plotAtTime(T_t, Nx, Ny, lx, ly, "t="+str(t_max)+" s")
