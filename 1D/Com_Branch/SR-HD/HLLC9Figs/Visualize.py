import numpy as np
import matplotlib.pyplot as plt

File1 = "0_0_Data"
File2 = "0_9_Data"
File3 = "0_99_Data"

File4 = "9_0_Data"
File5 = "9_9_Data"
File6 = "9_99_Data"

File7 = "99_0_Data"
File8 = "99_9_Data"
File9 = "99_99_Data"

x = np.arange(0,1,1/400)

Files = [[File1,File2,File3],[File4,File5,File6],[File7,File8,File9]]

fig, ax = plt.subplots(3,3, sharey=True, sharex=True)
LTitles = ["Vy_L=0.0,", "Vy_L=0.9,", "Vy_L=0.99,"]
RTitles = ["Vy_R=0.0", "Vy_R=0.9", "Vy_R=0.99"]

for i in range(3):
  for j in range(3):
    Dens = np.loadtxt(Files[i][j] + "/Density.dat")
    Xvel = np.loadtxt(Files[i][j] + "/VelocityX.dat")
    Yvel = np.loadtxt(Files[i][j] + "/VelocityY.dat")
    Pres = np.loadtxt(Files[i][j] + "/Pressure.dat")
    # rcm = np.loadtxt(Files[i][j] + "/Rcm.dat")

    ax[i][j].plot(x,Dens/25,'k-')
    ax[i][j].plot(x,Xvel,'b-')
    ax[i][j].plot(x,Yvel,'g-')
    ax[i][j].plot(x,Pres/1000,'r-')
    # ax[i][j].scatter(x,rcm*.5)
    ax[i][j].set_title("Weno, HLLC, RK3, " + LTitles[i] + " " + RTitles[j] )



plt.show()
