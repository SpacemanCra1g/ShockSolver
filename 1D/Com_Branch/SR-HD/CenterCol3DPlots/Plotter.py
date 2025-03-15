import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

FilePath = "./Pressure.dat"

var = np.loadtxt(FilePath)


plt.style.use('_mpl-gallery')
X = np.arange(0,1,1/8000)
Y = np.arange(0,1,1/30)

X,Y = np.meshgrid(X,Y)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(X, Y, var/max(var[0]), vmin=var.min() * 2, cmap=cm.jet)
ax.set_ylabel("V_y")
ax.set_xlabel("X")
ax.set_zlabel("Pressure")
plt.show()
#
# plt.plot(X,var[0])
# plt.show()





