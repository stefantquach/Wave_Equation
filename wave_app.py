import numpy as np
import tkinter as tk
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

def main():
    Lx=10
    Ly=10
    dx=0.1
    dy=0.1
    nx=int(Lx/dx)
    ny=int(Ly/dy)
    print(nx)
    print(ny)
    x=np.linspace(0,Lx,nx)
    y=np.linspace(0,Ly,ny)
    X, Y = np.meshgrid(x,y)

    T=10

    CFL = 0.5
    c = 1
    dt = CFL*dx/c

    K=0

    ## Actual simulation
    h = np.zeros([nx,ny])
    h_prev = h
    h_next = h

    t = 0
    fig = plt.figure()
    ax = mplot3d.Axes3D(fig)
    while(t < T):
        # Reflecting boundary condition
        h[:][[0, ny-1]]=0
        h[[0, nx-1]][:]=0

        # updating
        t += dt
        h_prev = h
        h = h_next

        if(t<10):
            h[49][49]=dt*dt*20*np.sin(30*np.pi*t/20)

        for i in range(1,nx-1):
            for j in range(1, ny-1):
                h_next[i][j] = 1/(1+K*dt)*(2*h[i][j] - h_prev[i][j] + K*dt*h_prev[i][j] + CFL*CFL *(h[i+1][j]+h[i][j+1]-4*h[i,j]+h[i-1][j]+h[i][j-1]));

        # ax = plt.axes(projection='3d')
        ax.cla()
        # ax.plot_surface(X,Y,h,cmap='viridis', edgecolor='none')
        # ax.scatter(X,Y,h)
        ax.set_xlim3d(0, Lx)
        ax.set_ylim3d(0, Ly)
        ax.set_zlim3d(-0.05, 0.05)

        plt.pause(0.001)

if __name__ == '__main__':
    main()
