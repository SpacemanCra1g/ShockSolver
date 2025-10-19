import numpy as np
import matplotlib.pyplot as plt 

p = np.loadtxt("./OutputData/Pressure.dat")
E = np.loadtxt("./OutputData/Energy.dat")
MomX = np.loadtxt("./OutputData/MomX.dat")
MomY = np.loadtxt("./OutputData/MomY.dat")
MomZ = np.loadtxt("./OutputData/MomZ.dat")
ConRho = np.loadtxt("./OutputData/ConDensity.dat")
rho = np.loadtxt("./OutputData/Density.dat")
u = np.loadtxt("./OutputData/VelocityX.dat")
v = np.loadtxt("./OutputData/VelocityY.dat")
rcm = np.loadtxt("./OutputData/Rcm.dat")
DivP = np.loadtxt("./OutputData/DivP.dat")
# w = np.loadtxt("./OutputData/VelocityZ.dat")
tempU = np.loadtxt("./Density.dat")

RCMrho = np.loadtxt("./OutputData/Density.dat")
RCMu = np.loadtxt("./OutputData/VelocityX.dat")

# Hp = np.loadtxt("./HLLC_Comp/Pressure.dat")
# Hrho = np.loadtxt("./HLLC_Comp/Density.dat")
# Hu = np.loadtxt("./HLLC_Comp/VelocityX.dat")

# Ep = np.loadtxt("./ExactSolution/Pressure.dat")
# Erho = np.loadtxt("./ExactSolution/Density.dat")
# Eu = np.loadtxt("./ExactSolution/VelocityX.dat")


# plotVar = rho

# if len(np.shape(p)) == 1 or True:
if  True:
    with  open("include/Parameters.h",'r') as f:
        param = f.readlines()
    N = 0
    xstart = 0
    xend = 0
    VL = 0
    VR = 0
    RS = 0
    Method = 0

    i = 0
    while not N or not xstart or not xend or not VL or not VR or not RS:
        if "#define NX " in param[i]:
            N = i
        if "#define X0 " in param[i]:
            xstart = i
        if "#define XN " in param[i]:
            xend = i
        if "#define YVELL " in param[i]:
            VL = i
        if "#define YVELR " in param[i]:
            VR = i
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
    VL = str(param[VL][14:-1])
    VR = str(param[VR][14:-1])
    RS = str(param[RS][16:-1])
    Method = str(param[Method][20:-1])

    if RS == "AUSMPLUS":
        RS = "AUSM+"
    elif RS == "AUSMPLUSUP":
        RS = "AUSM+up"


    Indx = np.argmax(p)
    # print(p[Indx-1:Indx+2])
    # print(E[Indx-1:Indx+2])
    # print(Indx)
    # exit()



    deltaX = (xend-xstart)/N
    xstart += deltaX/2 #adjust the interval one half deltax away from the start

    x = np.arange(xstart,xend,deltaX)

    # deltaX = (xend-xstart)/len(Erho)
    # Ex = np.arange(xstart,xend,deltaX)

    # deltaX = (xend-xstart)/len(Hrho)
    # Hx = np.arange(xstart,xend,deltaX)

    # print(x)



    DivP -= 200*deltaX
    DivP/=max(DivP)
    # DivP/=25
    # plt.plot(x,u-.3,'b-',linewidth=3)
    # plt.plot(x,v,'g-',linewidth=3)
    # plt.plot(x,rho/25,'k-',linewidth=3)
    # plt.plot(x,p/1000,'r--',linewidth=3)
    # plt.grid()
    # plt.plot(x,(DivP - 10*deltaX)/max(abs(DivP - 10*deltaX)) ,'y-')
    # plt.plot(x,DivP,'y-')
    # plt.scatter(x,rcm*.5)

    # plt.plot(Hx,Hu,'b.')
    # plt.plot(Hx,Hrho/25,'k.')
    # plt.plot(Hx,Hp/1000,'r.')

    # plt.plot(Ex,Eu,'b--')
    # plt.plot(Ex,Erho/25,'k--')
    # plt.plot(Ex,Ep/1000,'r--')
    # plt.plot(x,(w*w + u*u + v*v),'g')
    # plt.scatter(x,u,color='b',s=5, marker='.')
    # plt.scatter(x,rho/25,color='k',s=5,marker='.')
    # plt.scatter(x,p/1000,color='r',s=5,marker='.')
    title = "t = 0.4" +", V_yL = " + VL + ", V_yR = " + VR + ", Nx = " + str(N) + ", Contact Only"

    # plt.plot(x,rho,'k.-')
    # fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    # fig.suptitle("Contact wave only",fontsize=20 )
    # ax1.plot(x,rho/25,'k.-')
    # ax1.plot(x,u,'b.-')
    # ax2.plot(x,RCMrho/25,'k.-')
    # ax2.plot(x,RCMu,'b.-')
    # ax1.grid()
    # ax2.grid()
    # ax1.set_title("First Order Godunov",fontsize=16)
    # ax2.set_title("RCM",fontsize=16)
    # fig.legend(["Density","Normal Velocity"],fontsize=16)
    # plt.title("Shock tube problem at $t = 0.4$")
    # plt.legend(["Normal Velocity","Transverse Velocity","Density","Pressure"])

    # print(test.index(min(test)))
    # print(test)
    # print(p[199:202])
    fig, (ax2,ax3) = plt.subplots(1, 2, sharey=False)
    fig.suptitle("Contact wave only. 20 Timesteps in: FOG, HLLC",fontsize=20 )

    # ax1.plot(x,ConRho,'k',linewidth=3)
    # ax1.plot(x,E,'r',linewidth=3)
    # ax1.plot(x,MomX,'b',linewidth=3)
    # ax1.plot(x,MomY,'g--',linewidth=3)
    # ax1.legend(["Rho","Energy","MomentumX","MomentumY"])
    # ax1.set_title("Conservative Variables")
    # # ax1.set_xlim([.45,.55])
    # ax1.grid()

    # print(ConRho[200],MomX[200],MomY[200],E[200])

    # fig.set_xlim([.48,.52])

    ax2.plot(x,rho/25,'k.',linewidth=3)
    ax2.plot(x,p,'r^',linewidth=3,markersize=6,markerfacecolor='none')
    ax2.plot(x,u,'b.',linewidth=3)
    ax2.plot(x,v,'g.',linewidth=3)
    ax2.legend(["Rho/25","Pressure","Vx","Vy"])
    # ax2.set_xlim([.45,.55])
    ax2.grid()
    ax2.set_title("Primative Variables")
    # print(max(v))
    # print(v)

    gam = 5/3
    sig = gam/(gam-1)
    LorD = [1.0 / np.sqrt(1 - v[i]**2 - u[i]**2) for i in range(len(x))]
    hD = [1 + sig*(p[i]/rho[i]) for i in range(len(x))]
    DD = [LorD[i]*rho[i] for i in range(len(x))]
    SiD = [(LorD[i]**2)*rho[i]*hD[i]*u[i] for i in range(len(x))]
    SjD = [(LorD[i]**2)*rho[i]*hD[i]*v[i] for i in range(len(x))]
    ED = [(LorD[i]**2)*rho[i]*hD[i] - p[i] for i in range(len(x))]

    ax3.plot(x,DD,'k.',linewidth=3)
    ax3.plot(x,ED,'r^',linewidth=3,markersize=6,markerfacecolor='none')
    ax3.plot(x,SiD,'b.',linewidth=3)
    ax3.plot(x,SjD,'g.',linewidth=3)
    
    ax3.legend(["Rho","Energy","MomentumX","MomentumY"])
    ax3.set_title("Conservative Variables")
    # ax2.set_xlim([.45,.55])
    ax3.grid()
    

    plt.show()

    # plt.legend([RS +" Vx",RS +" Rho",RS+" Pres", "HLLC Vx","HLLC Rho","HLLC Pres","Exact Vx","Exact Rho","Exact Pres",])
    # plt.legend(["Vx",'Vy'," Rho","Pres","DivP"])
    # plt.legend(["Rho"])
    #
    # fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

# First plot
#     ax1.plot(x, tempU,"k.")
#     ax1.set_title("HLLC+WENO")
#     ax1.set_ylabel("Density",fontsize=14)
#     ax1.set_xlim([.45,.8])
#     ax1.grid()
#     ax1.legend()

# # Second plot
#     ax2.plot(x, rho, "k.")
#     ax2.set_title("RCM")
#     ax2.grid()
#     ax2.set_xlim([.45,.8])
#     ax2.legend()
#     fig.suptitle("$V_L=0.9 = V_R$ Shockwave only",fontsize=15)

# # Adjust layout
#     plt.tight_layout()
#     plt.show()

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


