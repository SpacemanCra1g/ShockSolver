#!/usr/bin/env python3

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.integrate as INT
import scipy.optimize as opt
r = 0
v = 1
vt = 2
p = 3
gamma = 5/3
RR = 0
RS = 1
SS = 2


def Geth(p,r):
    return 1 + (gamma/(gamma-1))*p/r

def Solve_RR_vx(state,pres,sign):
    # Compute Entropy
    SL = state[p]/(state[r]**gamma)
    A = Geth(state[p],state[r])*(1/np.sqrt(1 - (state[v]**2 + state[vt]**2) ))*state[vt]

    def integrate(p):
        rho = (p/SL)**(1/gamma)
        h = Geth(p,rho)
        cs = np.sqrt(gamma*p/(h*rho))
        
        return np.sqrt(h*h + A*A*(1-cs*cs))/((h*h + A*A)*rho*cs)
    

    # Rarefacetion Wave V^x_b Equation 4.201
    B1 = .5*np.log( (1 + state[v])/(1-state[v]))
    B2 = INT.quad(integrate,state[p],pres)[0]
    
    return np.tanh(B1+sign*B2)


def Solve_RR(Left,Right,p):
    ux3 = Solve_RR_vx(Left,p,-1)
    ux4 = Solve_RR_vx(Right,p,1)

    v13=GetRelSpeed(Left[v],ux3)
    v24=GetRelSpeed(Right[v],ux4)

    return GetRelSpeed(v13,v24)
    
        
def GetH(State):
    return 1 + gamma*State[p]/(State[r]* (gamma-1))

def GetA(State):
    return (1/np.sqrt(1-State[vt]**2 + State[v]**2))*State[vt]*GetH(State)

def GetRelSpeed(Left,Right): #4.101
    return (Left - Right)/(1 - Left*Right)

def GetVt(A,vx,h):
    return np.sqrt(A*A*(1- vx*vx)/(h*h + A*A))

def FindRiemannWave(State,ustar,pstar):
    S = State[p]/(State[r]**gamma)
    A = Geth(State[p],State[r])*((1 - (State[v]**2 + State[vt]**2))**(-.5))*State[vt]
    rho = (pstar / S)**(1/gamma)
    vT = GetVt(A,ustar,Geth(pstar,rho))
    return [rho,ustar,vT,pstar]
    
    

def RightWaveType(StateL,StateR):
    v12_0 = GetRelSpeed(StateL[v],StateR[v])
    A1 = GetA(StateL)
    S1 = StateL[p]/(StateL[r]**(gamma))
    rFp = lambda p: (p/S1)**(1/gamma)
    hFp = lambda p,rho: Geth(p,rho)
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
    hFp2 = lambda p,rho: Geth(p,rho)
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
        eps = 1e-15
        p_min = (StateR[p] + eps)*eps
        p_max = StateL[p]
        p_star = opt.brentq(lambda p: Solve_RR(StateL,StateR,p) - v12_0, p_min,p_max)
        print(p_star)
        vstar= Solve_RR_vx(StateR,p_star,1)
        print(vstar)

        Wave3 = FindRiemannWave(StateL,vstar,p_star)
        Wave3Tick = FindRiemannWave(StateR,vstar,p_star)
        print(Wave3)
        print(Wave3Tick)
        
        

    else:
       h1 = GetH(StateL[p],StateL[r])
       h2 = GetH(StateR[p],StateR[r])
       
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
    StateL = [1, 0, .9, 1000]
    StateR = [1, 0,.9,.01]
    RightWaveType(StateL,StateR)
    
    
