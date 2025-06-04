#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

p = np.loadtxt("./Pres.dat")
x = np.linspace(1.e-3,200,len(p))


plt.plot(x,p)
plt.show()
