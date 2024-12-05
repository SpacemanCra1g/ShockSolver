import numpy as np
import matplotlib.pyplot as plt 

p = np.loadtxt("./OutputData/Pressure.dat")
rho = np.loadtxt("./OutputData/Density.dat")
u = np.loadtxt("./OutputData/VelocityX.dat")

plotVar = rho

if len(np.shape(p)) == 1:
    with  open("include/Parameters.h",'r') as f:
        param = f.readlines()
    N = 0
    xstart = 0
    xend = 0
    RS = 0
    i = 0
    Method = 0
    while not N or not xstart or not xend or not RS:
        if "#define NX " in param[i]:
            N = i
        if "#define X0 " in param[i]:
            xstart = i
        if "#define XN " in param[i]:
            xend = i
        if "#define RIEMANN " in param[i]:
            RS = i
        if "#define SpaceMethod " in param[i]:
            Method = i
        i += 1
        if i == len(param):
            print("Coefficient not found in file\n Exiting")
            exit(0)

    N = int(param[N][11:-1])
    xstart = float(param[xstart][11:-1])
    xend = float(param[xend][11:-1])
    RS = str(param[RS][16:-1])
    Method = str(param[Method][20:-1])


    deltaX = (xend-xstart)/N
    xstart += deltaX/2 #adjust the interval one half deltax away from the start

    x = np.arange(xstart,xend,deltaX)
    # print(x)

    maxu = max(u);
    maxp = max(p);
    # plt.plot(x,u/maxu,'b')
    plt.plot(x,rho,'k-')
    # plt.plot(x,p/maxp,'r')
    # plt.plot(x,(w*w + u*u + v*v),'g')
    # plt.scatter(x,u,color='b',s=5, marker='.')
    # plt.scatter(x,rho/25,color='k',s=5,marker='.')
    # plt.scatter(x,p/1000,color='r',s=5,marker='.')
    # title =  "SpaceMethod = " + Method +", Nx = " + str(N) + ", " + RS
    title = "U2 Turned off, ShuOsher"
    plt.title(title)
    plt.grid()
    # plt.legend(["Vx","Rho","Pres"])

    plt.show()

elif len(np.shape(rho)) == 2:
    with  open("Parameters",'r') as f:
        param = f.readlines()

    Nx = int(param[13][2:])
    xstart = float(param[14][2:])
    xend = float(param[15][2:])

    Ny = int(param[20][2:])
    ystart = float(param[21][2:])
    yend = float(param[22][2:])

    deltaX = (xend-xstart)/Nx
    xstart += deltaX/2 #adjust the interval one half deltax away from the start

    x = np.arange(xstart,xend,deltaX)

    deltaY = (yend-ystart)/Ny
    xstart += deltaY/2 #adjust the interval one half deltax away from the start

    y = np.arange(ystart,yend,deltaY)

    X,Y = np.meshgrid(x,y)


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.pcolormesh(X,Y,rho, cmap='seismic')
    fig.colorbar(ax.pcolormesh(X, Y, plotVar, cmap='seismic'), orientation="horizontal" )
    # ax.axis([x.min(),x.max(),y.min(),y.max()])

    # fig,ax = plt.subplots(subplot_kw={"projection":"3d"})
    # c = ax.plot_surface(X,Y,rho)

    ax.set_title('NN Order 5')

    plt.show()


