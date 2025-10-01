import numpy as np
import matplotlib.pyplot as plt 

p = np.loadtxt("./OutputData/Pressure.dat")
rho = np.loadtxt("./OutputData/Density.dat")
u = np.loadtxt("./OutputData/VelocityX.dat")

X = [[0.5]*4,[0.204196, 0.482432, 0.731863, 0.938039]]
Y = [[0]*4,[.25]*4]

N = len(p)
xend = 1
xstart = 0
deltaX = (xend-xstart)/N
xstart += deltaX/2 #adjust the interval one half deltax away from the start

x = np.arange(xstart,xend,deltaX)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True)
# plt.plot(x,u,'b-')
ax1.plot(x,rho,'k-')
colors = ["r","r","k","b"]
labels = ["$\mathcal{W}_1^h$","$\mathcal{W}_1^t$","$\mathcal{C}$","$\mathcal{W}_2$"]
for i in range(4):
    ax1.plot([X[1][i],X[1][i]], [Y[0][i],Y[1][i] + .75],color=colors[i])
# plt.plot(x,p,'r-')



for i in range(4):
    ax2.plot([X[0][i],X[1][i]], [Y[0][i],Y[1][i]],color=colors[i],label=labels[i])
ax2.set_xlim(0,1)
ax2.set_ylim(0,.25)
ax2.grid()
ax1.grid()
ax2.legend()
ax1.set_ylim(0,1.05)
ax2.set_ylabel("Time")
ax2.set_xlabel("Domain")
ax1.set_ylabel("Density")
fig.suptitle("Sod shock tube at time $t = 0.25$")
ax2.set_title("Sod shock tube, Riemann Fan")
plt.show()
