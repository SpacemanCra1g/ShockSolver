import numpy as np
import matplotlib.pyplot as plt 

p = np.loadtxt("./OutputData/Pres.dat")
rho = np.loadtxt("./OutputData/Density.dat")
u = np.loadtxt("./OutputData/XVel.dat")
v = np.loadtxt("./OutputData/YVel.dat")


plotVar = rho[3,3:-3]

if len(np.shape(p)) == 1 or 1:
    with  open("Parameters",'r') as f:
        param = f.readlines()

    N = int(param[14][2:])
    xstart = float(param[15][2:])
    xend = float(param[16][2:])

    deltaX = (xend-xstart)/N
    xstart += deltaX/2 #adjust the interval one half deltax away from the start

    x = np.arange(xstart,xend,deltaX)

    plt.plot(x,plotVar,'r-', markerfacecolor='none')
    plt.title("GP-R1 RK3, nx = 256 NN")
    plt.grid()

    plt.show()

elif len(np.shape(rho)) == 2:
    with  open("Parameters",'r') as f:
        param = f.readlines()

    Nx = int(param[14][2:])
    xstart = float(param[15][2:])
    xend = float(param[16][2:])

    Ny = int(param[21][2:])
    ystart = float(param[22][2:])
    yend = float(param[23][2:])

    deltaX = (xend-xstart)/Nx
    xstart += deltaX/2 #adjust the interval one half deltax away from the start

    x = np.arange(xstart,xend,deltaX)

    deltaY = (yend-ystart)/Ny
    xstart += deltaY/2 #adjust the interval one half deltax away from the start

    y = np.arange(ystart,yend,deltaY)

    X,Y = np.meshgrid(x,y)


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    # ax.pcolormesh(X,Y,rho, cmap='seismic')
    # fig.colorbar(ax.pcolormesh(X, Y, plotVar, cmap='seismic'), orientation="horizontal" )
    ax.axis([x.min(),x.max(),y.min(),y.max()])

    fig,ax = plt.subplots(subplot_kw={"projection":"3d"})
    c = ax.plot_surface(X,Y,rho)

    ax.set_title('NN Order 5')

    plt.show()


