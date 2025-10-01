import numpy as np
import matplotlib.pyplot as plt 


HLLCrho = np.loadtxt("./HLLCDensity.dat")
RCMrho = np.loadtxt("./RCMDensity.dat")

OverlayExact = False
PathToExact = "./ExactSolutions/SlowShockTORO/"


if OverlayExact:
    pE = np.loadtxt(PathToExact + "Pressure.dat")
    rhoE = np.loadtxt(PathToExact + "Density.dat")
    uE = np.loadtxt(PathToExact + "VelocityX.dat")



if True:
    with  open("include/Parameters.h",'r') as f:
        param = f.readlines()
    N = 0
    xstart = 0
    xend = 0
    RS = 0
    i = 0
    Method = 0
    Problem = 0
    while not N or not xstart or not xend or not RS or not Problem:
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
        if "#define TestProblem " in param[i]:
            Problem = i
        i += 1
        if i == len(param):
            print("Coefficient not found in file\n Exiting")
            exit(0)

    N = int(param[N][11:-1])
    xstart = float(param[xstart][11:-1])
    xend = float(param[xend][11:-1])
    RS = str(param[RS][16:-1])
    if RS == "AUSMPLUS":
        RS = "AUSM+"
    if RS == "AUSMPLUSUP":
        RS = "AUSM+up"
    Method = str(param[Method][20:-1])
    Problem = str(param[Problem][20:-1])



    deltaX = (xend-xstart)/N
    xstart += deltaX/2 #adjust the interval one half deltax away from the start

    x = np.arange(xstart,xend,deltaX)

    if OverlayExact:
        xE = np.linspace(xstart,xend,num=len(rhoE),endpoint=True)
    # print(x)

    
    # plt.plot(x,rho,'ko',markersize=2)

    fig, (ax1,ax2) = plt.subplots(1,2,sharey=True)
    ax1.plot(x,HLLCrho,'k.')
    ax2.plot(x,RCMrho,'k.')
    
    
    # fig.legend(["Density"])
    
    title =  "Hydro slow shock at T = 1.0"
    ax1.set_title("HLLC+WENO")
    ax2.set_title("RCM")
    ax1.set_xlim([-6,2.5])
    ax2.set_xlim([-6,2.5])
    # title = "U2 Turned off, ShuOsher"
    fig.suptitle(title,fontsize=14)
    ax1.grid()
    ax2.grid()
    # plt.legend(["Rho","Vx","Pres"])
    fig.show()
    fig.tight_layout()
    fig.legend(["Density"],loc="center right")
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


