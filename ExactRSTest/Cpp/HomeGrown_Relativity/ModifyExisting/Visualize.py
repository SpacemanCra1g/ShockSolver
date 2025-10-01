import numpy as np
import matplotlib.pyplot as plt 

p = np.loadtxt("./OutputData/Pressure.dat")
rho = np.loadtxt("./OutputData/Density.dat")
u = np.loadtxt("./OutputData/VelocityX.dat")
v = np.loadtxt("./OutputData/VelocityY.dat")

if max(p) > 100:
    p /= 1000
    rho /= 25


N = len(p)
xend = 1
xstart = 0
deltaX = (xend-xstart)/N
xstart += deltaX/2 #adjust the interval one half deltax away from the start

x = np.arange(xstart,xend,deltaX)


# plt.plot(x,u,'b-')
# plt.plot(x,v,'g-')
plt.plot(x,rho,'k-')
# plt.plot(x,p,'r-')

plt.grid()
# plt.legend(["Vx","Vy"," Rho","Pres"])
# plt.legend(["Rho"])

plt.show()


