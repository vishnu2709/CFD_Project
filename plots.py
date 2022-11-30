import matplotlib.pyplot as plt
import numpy as np

pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'
pl.rcParams['xtick.major.size'] = 8
pl.rcParams['xtick.minor.size'] = 4
pl.rcParams['xtick.major.pad']  = 8
pl.rcParams['xtick.minor.pad']  = 8
pl.rcParams['xtick.color']      = 'k'
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'
pl.rcParams['ytick.major.size'] = 8
pl.rcParams['ytick.minor.size'] = 4
pl.rcParams['ytick.major.pad']  = 8
pl.rcParams['ytick.minor.pad']  = 8
pl.rcParams['ytick.color']      = 'k'
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'

# You may add additional functions for ploting as you see fit

def VelocityFiled(x,y,u,v):
    xx, yy = np.meshgrid(x,y)
    plt.streamplot(xx,yy,u, v, color=np.sqrt(u*u + v*v),density=1.5,linewidth=1.5, cmap=plt.cm.viridis)
    plt.colorbar(label = 'velocity [m/s]')
    plt.xlabel('x [m]',fontsize = 14 )
    plt.ylabel('y [m]',fontsize = 14 )
    plt.title('Streamlines', fontsize = 14)
    plt.tick_params(labelsize=12)
    plt.ylim([0,0.04])
    plt.xlim([0,0.04])
    plt.show()


def plotAtTime(T, Nx, Ny, Lx, Ly, plot_title):
    x = np.linspace(0, Lx, Nx+1)
    y = np.linspace(0, Ly, Ny+1)
    X, Y = np.meshgrid(x, y)
    pl.pcolormesh(X, Y, T)
    pl.title(plot_title)
    pl.xlabel("$x$ (ft)")
    pl.ylabel("$y$ (ft)")
    pl.colorbar(label="Temperature (K)")
    pl.savefig(plot_title+".pdf")
    pl.show()

def plotLines(T, N, L, times, plot_title):
    x = np.linspace(0, L, N+1)
    for i in range(len(T)):
        pl.plot(x, T[i], label="t="+str(times[i])+" hr")
    pl.xlabel("$y (ft)$")
    pl.ylabel("$T (K)$")
    pl.legend()
    pl.title(plot_title)
    pl.savefig("line-comparison.pdf")
    pl.show()
