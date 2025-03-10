#!/usr/bin/env python3

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.integrate as INT
r = 0
v = 1
vt = 2
p = 3
gamma = 5/3
RR = 0
RS = 1
SS = 2


def GetH(State):
    return 1 + gamma*State[p]/State[r]

def GetA(State):
    return (1/np.sqrt(1-State[vt]**2 + State[v]**2))*State[vt]*GetH(State)

def GetRelSpeed(Left,Right): #4.101
    return (Left - Right)/(1 - Left*Right)

def RightWaveType(StateL,StateR):
    v12_0 = GetRelSpeed(StateL[v],StateR[v])
    A1 = GetA(StateL)
    S1 = StateL[p]/(StateL[r]**(gamma))
    rFp = lambda p: (p/S1)**(1/gamma)
    hFp = lambda p,rho: 1 + gamma*p/rho
    csFp = lambda p,h,rho: np.sqrt(gamma*p/(h*rho))
    lor1 = (1 - (StateL[v]**2 + StateL[vt]**2))**(-.5)
    def Int1(p):
        r = rFp(p)
        h = hFp(p,r)
        cs = csFp(p,h,r)
        # A1 = h*lor1*StateL[vt]
        return np.sqrt(h**2 + (A1**2)*(1 - cs**2))/( (h**2 + A1**2)*r*cs)
    v1_x = np.tanh(INT.quad(Int1,StateL[p],0)[0]) # Equation 4.214
    

    A2 = GetA(StateR)
    S2 = StateR[p]/(StateR[r]**(gamma))
    rFp2 = lambda p: (p/S2)**(1/gamma)
    hFp2 = lambda p,rho: 1 + gamma*p/rho
    csFp2 = lambda p,h,rho: np.sqrt(gamma*p/(h*rho))
    lor2 = (1 - (StateR[v]**2 + StateR[vt]**2))**(-.5)
    def Int2(p):
        r = rFp2(p)
        h = hFp2(p,r)
        cs = csFp2(p,h,r)
        # A2 = h*lor2*StateR[vt]
        return np.sqrt((h**2 + (A2**2)*(1 - cs**2)))/( (h**2 + A2**2)*r*cs)
    
    v2_x = np.tanh(INT.quad(Int2,0,StateR[p])[0]) # Equation 4.215
    # print(v2_x)
    # print(v1_x)
    # exit()
    v_12x_2R = GetRelSpeed(v1_x,v2_x)
    v_12x_SR = np.tanh(INT.quad(Int1,StateL[p],StateR[p])[0])
    
    if v12_0 <= v_12x_SR:
        print("Two Rarefaction Case")

    else:
       h1 = 1 + gamma*StateL[p]/StateL[r]
       h2 = 1 + gamma*StateR[p]/StateR[r]
       
       D =  4*gamma*StateL[p]*( ( (gamma-1)*StateR[p] + StateL[p])/ ( ((gamma-1)*(StateL[p] - StateR[p]))**2   ))  
       D*= (h2*(StateR[p] - StateL[p])/StateR[r] - h2*h2)
       D = 1 - D
       
       h3 = (np.sqrt(D) - 1)*(gamma-1)*(StateL[p] - StateR[p])/( 2*( (gamma-1)*StateR[p] + StateL[p]))
       
       J2 = -(gamma/(gamma-1))*(StateL[p] -StateR[p])
       J2 /= ( h3*(h3-1)/StateL[p] - h2*(h2-1)/StateR[p])

       lor2 = 1/(1 - StateR[v]**2 - StateR[vt]**2)
       Vs = ((StateR[r]**2)*lor2*StateR[v] + np.sqrt(np.abs(J2))*np.sqrt(J2 + (StateR[r]**2)*lor2*(1 - StateR[v]**2)))
       Vs /= ((StateR[r]**2)*lor2 + J2)
       
       v_12x_SS = (StateL[p] - StateR[p])*(1- StateR[v]*Vs)
       v_12x_SS /= (Vs- StateR[v])*(h2*StateR[r]*lor2*(1-StateR[v]**2) + StateL[p] - StateR[p])


       if v12_0 <= v_12x_SS:
           print("One Shock and one Rarefaction")
           print(Vs)
           print(v12_0)
       else:
           print("Two Shocks")
           print(Vs)
           print(v_12x_SS)
           print(v12_0)
       
if __name__ == "__main__":
    StateL = [1, 0, .999, 1]
    StateR = [.125, .5,.0,.1]
    RightWaveType(StateL,StateR)
        
    
