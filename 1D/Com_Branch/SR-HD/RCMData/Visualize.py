import numpy as np
import matplotlib.pyplot as plt 

p = np.loadtxt("./0_0Data/Pressure.dat")
rho = np.loadtxt("./0_0Data/Density.dat")
u = np.loadtxt("./0_0Data/VelocityX.dat")
v = np.loadtxt("./0_0Data/VelocityY.dat")

xend = 1.0
xstart = 0.0
N = 400
deltaX = (xend-xstart)/N
xstart += deltaX/2

x = np.arange(xstart,xend,deltaX)
fig,ax = plt.subplots(3,3, sharex = True, sharey= True)
# ax1 = plt.subplot(911)
# ax2 = plt.subplot(912)
# ax3 = plt.subplot(913)
# ax4 = plt.subplot(921)
# ax5 = plt.subplot(922)
# ax6 = plt.subplot(923)
# ax7 = plt.subplot(931)
# ax8 = plt.subplot(932)
# ax9 = plt.subplot(933)
# ax = [[ax1,ax2,ax3],[ax4,ax5,ax6],[ax7,ax8,ax9]]


# Figure 1
ax[0,0].plot(x,u,'b-')
ax[0,0].plot(x,v,'g-')
ax[0,0].plot(x,rho/25,'k-')
ax[0,0].plot(x,p/1000,'r-')

title = "V_yL = 0.0, V_yR = 0.0"
ax[0,0].set_title(title)

# Figure 2
p = np.loadtxt("./0_9Data/Pressure.dat")
rho = np.loadtxt("./0_9Data/Density.dat")
u = np.loadtxt("./0_9Data/VelocityX.dat")
v = np.loadtxt("./0_9Data/VelocityY.dat")
title = "V_yL = 0.0, V_yR = 0.9"

ax[0,1].plot(x,u,'b-')
ax[0,1].plot(x,v,'g-')
ax[0,1].plot(x,rho/25,'k-')
ax[0,1].plot(x,p/1000,'r-')
ax[0,1].set_title(title)

# Figure 3
p = np.loadtxt("./0_99Data/Pressure.dat")
rho = np.loadtxt("./0_99Data/Density.dat")
u = np.loadtxt("./0_99Data/VelocityX.dat")
v = np.loadtxt("./0_99Data/VelocityY.dat")
title = "V_yL = 0.0, V_yR = 0.99"

ax[0,2].plot(x,u,'b-')
ax[0,2].plot(x,v,'g-')
ax[0,2].plot(x,rho/25,'k-')
ax[0,2].plot(x,p/1000,'r-')
ax[0,2].set_title(title)

# Figure 4
p = np.loadtxt("./9_0Data/Pressure.dat")
rho = np.loadtxt("./9_0Data/Density.dat")
u = np.loadtxt("./9_0Data/VelocityX.dat")
v = np.loadtxt("./9_0Data/VelocityY.dat")
title = "V_yL = 0.9, V_yR = 0.0"

ax[1,0].plot(x,u,'b-')
ax[1,0].plot(x,v,'g-')
ax[1,0].plot(x,rho/25,'k-')
ax[1,0].plot(x,p/1000,'r-')
ax[1,0].set_title(title)

# Figure 5
p = np.loadtxt("./9_9Data/Pressure.dat")
rho = np.loadtxt("./9_9Data/Density.dat")
u = np.loadtxt("./9_9Data/VelocityX.dat")
v = np.loadtxt("./9_9Data/VelocityY.dat")
title = "V_yL = 0.9, V_yR = 0.9"

ax[1,1].plot(x,u,'b-')
ax[1,1].plot(x,v,'g-')
ax[1,1].plot(x,rho/25,'k-')
ax[1,1].plot(x,p/1000,'r-')
ax[1,1].set_title(title)

# Figure 6
p = np.loadtxt("./9_99Data/Pressure.dat")
rho = np.loadtxt("./9_99Data/Density.dat")
u = np.loadtxt("./9_99Data/VelocityX.dat")
v = np.loadtxt("./9_99Data/VelocityY.dat")
title = "V_yL = 0.9, V_yR = 0.99"

ax[1,2].plot(x,u,'b-')
ax[1,2].plot(x,v,'g-')
ax[1,2].plot(x,rho/25,'k-')
ax[1,2].plot(x,p/1000,'r-')
ax[1,2].set_title(title)

# Figure 7
p = np.loadtxt("./99_0Data/Pressure.dat")
rho = np.loadtxt("./99_0Data/Density.dat")
u = np.loadtxt("./99_0Data/VelocityX.dat")
v = np.loadtxt("./99_0Data/VelocityY.dat")
title = "V_yL = 0.99, V_yR = 0.0"

ax[2,0].plot(x,u,'b-')
ax[2,0].plot(x,v,'g-')
ax[2,0].plot(x,rho/25,'k-')
ax[2,0].plot(x,p/1000,'r-')
ax[2,0].set_title(title)

# Figure 8
p = np.loadtxt("./99_9Data/Pressure.dat")
rho = np.loadtxt("./99_9Data/Density.dat")
u = np.loadtxt("./99_9Data/VelocityX.dat")
v = np.loadtxt("./99_9Data/VelocityY.dat")
title = "V_yL = 0.99, V_yR = 0.9"

ax[2,1].plot(x,u,'b-')
ax[2,1].plot(x,v,'g-')
ax[2,1].plot(x,rho/25,'k-')
ax[2,1].plot(x,p/1000,'r-')
ax[2,1].set_title(title)

# Figure 9
p = np.loadtxt("./99_99Data/Pressure.dat")
rho = np.loadtxt("./99_99Data/Density.dat")
u = np.loadtxt("./99_99Data/VelocityX.dat")
v = np.loadtxt("./99_99Data/VelocityY.dat")
title = "V_yL = 0.99, V_yR = 0.99"

ax[2,2].plot(x,u,'b-')
ax[2,2].plot(x,v,'g-')
ax[2,2].plot(x,rho/25,'k-')
ax[2,2].plot(x,p/1000,'r-')
ax[2,2].set_title(title)

for i in range(3):
    for j in range(3):
        ax[i,j].grid()


ax[0,1].sharey(ax[0,0])
ax[0,2].sharey(ax[0,0])


ax[0,0].legend(["Vx","Vy","Rho","Pres"])
fig.suptitle("T = 0.4, RK1, FOG, RCM",fontsize=18)

plt.show()
