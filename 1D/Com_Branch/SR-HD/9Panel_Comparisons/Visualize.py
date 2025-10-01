import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
Dir = ["ExactData/","HybridData/","HLLCData/"]
Lines=['','','-^']
Labels = ["Exact Solve","Hybrid Solve","HLLC Solve"]
line = [1,2,3]

for t in [0,2]:
  for i in range(3):
    for j in range(3,):
      File = Dir[t]+Files[i][j]
      Dens = np.loadtxt(File + "/Density.dat")
      Xvel = np.loadtxt(File + "/VelocityX.dat")
      Yvel = np.loadtxt(File + "/VelocityY.dat")
      Pres = np.loadtxt(File + "/Pressure.dat")
      # rcm = np.loadtxt(Files[i][j] + "/Rcm.dat")
      ax[i][j].grid()
      ax[i][j].plot(x,Pres/1000,'r-'+Lines[t],linewidth=2,markerfacecolor='none')
      # line[t], = ax[i][j].plot(x,Dens/25,Lines[t],linewidth=3,label=Labels[t])
      ax[i][j].plot(x,Xvel,'b-'+Lines[t],linewidth=2,markerfacecolor='none')
      ax[i][j].plot(x,Dens/25,'k-'+Lines[t],linewidth=2,markerfacecolor='none')
      # ax[i][j].plot(x,Yvel,'g')

      # ax[i][j].scatter(x,rcm*.5)
      ax[i][j].set_title("" + LTitles[i] + " " + RTitles[j] )
      ax[i][j].grid()

# print(line)
# plt.suptitle("Comparison of the exact solution v.s. HLLC+WENO",fontsize=20,y=.95)
# plt.legend(handles=[t for t in [0,1]],fontsize=18)

plt.show()


# Dens = np.loadtxt(Dir[0]+File8 + "/Density.dat")
# Vx = np.loadtxt(Dir[0]+File8 + "/VelocityX.dat")
# p = np.loadtxt(Dir[0]+File8 + "/Pressure.dat")
# plt.plot(x,Dens/25,'b--^')
# plt.plot(x,Vx,'b--^')
# plt.plot(x,p/1000,'b--^')

# Dens = np.loadtxt(Dir[1]+File8 + "/Density.dat")
# Vx = np.loadtxt(Dir[1]+File8 + "/VelocityX.dat")
# p = np.loadtxt(Dir[1]+File8 + "/Pressure.dat")
# plt.plot(x,Dens/25,'r-^')
# plt.plot(x,Vx,'r-')
# plt.plot(x,p/1000,'r-')
# # plt.xlim([.35,.65])



# color_legend = [
#     Line2D([0], [0], color='blue', lw=2, label='Exact Solution'),
#     Line2D([0], [0], color='red', lw=2, label='Hybrid Solution')
# ]

# marker_legend = [
#   Line2D([0], [0], color='black', marker='^', lw=0, label='Density'),
#     Line2D([0], [0], color='black', lw=2, label='X-Velocity')
# ]

# first_legend = plt.legend(handles=color_legend, title="Color meaning", loc='upper left')
# plt.gca().add_artist(first_legend)  # keep it when adding another legend
# plt.legend(handles=marker_legend, title="Marker meaning", loc='upper right')

# plt.title("V_L = 0.99, V_R = 0.9, Hybrid vs Exact. Density and $v^x$ ")
# plt.show()
