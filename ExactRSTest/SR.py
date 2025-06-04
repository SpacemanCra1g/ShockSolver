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
sigma = gamma/(gamma-1)
nnRR = 0
RS = 1
SS = 2

def SoundSpeed(State):
    h = 1 + sigma*State[p]/State[r]
    return np.sqrt(gamma*State[p]/(h*State[r]))

def Find_SS_Wave(StateL,StateR,vstar, pstar):
    Wave3 = FindShockWave(StateL,vstar,pstar)
    Wave3Tick = FindShockWave(StateR,vstar,pstar)
    return Wave3, Wave3Tick

def Find_RS_Wave(StateL, StateR,vstar,pstar):
    Wave3 = FindRiemannWave(StateL,vstar,pstar)
    Wave3Tick = FindShockWave(StateR,vstar,pstar)

    return Wave3, Wave3Tick


def FindShockWave(state,vstar,pstar):
    hB = Taub(state,pstar)
    A = Geth(state[p],state[r])*(1/np.sqrt(1 - (state[v]**2 + state[vt]**2) ))*state[vt]
    Wave = state.copy()
    Wave[p] = pstar
    Wave[v] = vstar
    Wave[vt] = A*np.sqrt( (1 - vstar**2)/(hB**2 + A**2) )
    Wave[r] = sigma/(hB-1)*pstar
    return Wave

def J_sqr(pres1,pres2,hA,hB):
    val = - sigma * (pres1 - pres2) / (
                hA * (hA - 1.) / pres1 - hB * (hB - 1.) / pres2)
    return val

def ShockSpeed(state,J,sign):
    lor = 1/np.sqrt(1 - state[v]**2 - state[vt]**2)
    D = state[r]*lor
    return (D ** 2 * state[v] + sign * J * np.sqrt(J ** 2 + D ** 2 * (1 - state[v] ** 2))) / (D ** 2 + J ** 2)

def Taub(state,pres):
    hA = Geth(state[p],state[r])
    c_2 = (1 + (state[p] - pres)/(pres*sigma) )
    c_1 = - (state[p] - pres)/(pres*sigma)
    c_0 = hA * (state[p] - pres)/state[r] - hA*hA

    if c_2 == 0:
        return -c_0/c_1
    else:
        return (-c_1 + np.sqrt(c_1*c_1 - 4 * c_2*c_0))/(2*c_2)

def Geth(p,r):
    return 1 + sigma*p/r

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

def Solve_Shock_vx(state,pres,sign):
    # Get H
    hA = Geth(state[p],state[r])
    # Taub Adiabat
    hB = Taub(state,pres)
    # print(hA, hB, hA-hB)
    J2 = J_sqr(state[p],pres,hA,hB)
    J = np.sqrt(np.abs(J2))
    Vs = ShockSpeed(state,J,sign)
    Ws = 1/(np.sqrt(1 - Vs*Vs) )
    lor = 1 / np.sqrt(1 - (state[v]**2 + state[vt]**2) )

    Value = (hA * lor * state[v] + sign * Ws * (pres - state[p]) / J) / (
                hA * lor + (pres - state[p]) * (
                    sign * Ws * state[v] / J + 1 / (state[r] * lor) ))
    # print("Shock Vel")
    # print(Value)
    return Value

def Solve_RR(Left,Right,p):
    ux3 = Solve_RR_vx(Left,p,-1)
    ux4 = Solve_RR_vx(Right,p,1)

    v13=GetRelSpeed(Left[v],ux3)
    v24=GetRelSpeed(Right[v],ux4)

    return GetRelSpeed(v13,v24)

def Solve_RS(StateL, StateR, p):
    ux3 = Solve_RR_vx(StateL,p,-1)
    ux4 = Solve_Shock_vx(StateR,p,1)
    # print(ux4)
    v13 = GetRelSpeed(StateL[v],ux3)
    v64 = GetRelSpeed(StateR[v],ux4)

    return GetRelSpeed(v13,v64)

def Solve_SS(StateL,StateR,p):
    ux3 = Solve_Shock_vx(StateL,p,-1)
    ux4 = Solve_Shock_vx(StateR,p,1)
    v13 = GetRelSpeed(StateL[v],ux3)
    v64 = GetRelSpeed(StateR[v],ux4)

    return GetRelSpeed(v13,v64)

        
def GetH(State):
    return 1 + gamma*State[p]/(State[r]* (gamma-1))

def GetA(State):
    return (1/np.sqrt(1-State[vt]**2 - State[v]**2)) *State[vt]*GetH(State)

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
    # v1_x = np.tanh(INT.quad(Int1,StateL[p],0)[0]) # Equation 4.214
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
        return Wave3, Wave3Tick
        
        

    else:
       h1 = GetH(StateL)
       h2 = GetH(StateR)
       
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

       # print (v12_0)
       # print (v_12x_SS)
       # exit()

       if v12_0 <= v_12x_SS:
           print("One Shock and one Rarefaction")

           eps = 1e-15
           p_min = StateR[p] + eps
           p_max = StateL[p]
           # print(Solve_RS(StateL, StateR, p_min) - v12_0)
           # print(Solve_RS(StateL, StateR, p_max) - v12_0)
           # exit(0)
           assert (p_min < p_max)
           p_star = opt.brentq(lambda p: Solve_RS(StateL, StateR, p) - v12_0, p_min, p_max)
           vstar = Solve_Shock_vx(StateR, p_star, 1)

           Wave3, Wave3Tick = Find_RS_Wave(StateL, StateR, vstar, p_star)

           # print(Wave3)
           # print(Wave3Tick)
           # print(Wave3[v]*.4)
           return Wave3, Wave3Tick

       else:
           print("Two Shocks")
           pstar = 1.1*StateL[p]
           pstar = opt.root(lambda p: Solve_SS(StateL,StateR,p) - v12_0, pstar).x[0]
           vstar = Solve_Shock_vx(StateR, pstar, 1)

           Wave3, Wave3Tick = Find_SS_Wave(StateL,StateR,vstar, pstar)

           print(Wave3)
           print(Wave3Tick)
           return Wave3, Wave3Tick


def RightShockVT(StateL,StateR,sign):
    state = StateL
    hA = Geth(state[p],state[r])
    # Taub Adiabat
    pres = StateR[p] #Wave3Tick[p]
    hB = Taub(state,pres)
    print("HB == ", hA)
    J2 = J_sqr(state[p],pres,hA,hB)
    J = np.sqrt(np.abs(J2))
    Vs = ShockSpeed(state,J,sign)
    print("Shock Location")
    print(Vs*.4 + .5)

def CsFromE(press,S):
   r = (press/S)**(1/gamma)
   h = 1 + sigma*press/r
   return np.sqrt(gamma*press/(h*r))

def SVel(state,S,sign):
    SP_sqr = state[v]**2 + state[vt]**2
    cs = CsFromE(state[p],S)
    cs_sqr = cs*cs
    return (state[v] * (1 - cs_sqr)
                + sign * cs * np.sqrt((1 - SP_sqr) * (1 - SP_sqr * cs_sqr - state[v] ** 2 * (1 - cs_sqr)))
                ) / (1 - SP_sqr * cs_sqr)

def RareFactionTails(StateL, StateR,sign):
    S = StateL[p]/(StateL[r]**gamma)
    S2 = StateR[p]/(StateR[r]**gamma)
    assert ( np.abs(S - S2) < 1e-10 *S )
    head = SVel(StateL,S,sign)
    tail = SVel(StateR,S,sign)
    print("Head velocity")
    print(head*.4 + .5)
    print("Tail velocity")
    print(tail*.4 + .5)

def ContactLocation(State):
    print("Contact Wave Location")
    print(State[v]*.4 + .5)

def ux(xi,StateL,pressure, sign, A):
    S = StateL[p]/np.power(StateL[r],gamma)
    rho = np.power(pressure/S, 1/gamma)
    h = 1 + sigma*pressure/rho
    cs = np.sqrt(gamma*pressure/(h*rho))
    a = cs*h
    b = sign*np.sqrt(A*A *(1 - cs*cs) + h*h)
    return (a-b*xi)/(a*xi - b)


def RarefactionState(xi,StateL,StateR,sign):
    pmin = min(StateL[p],StateR[p])
    pmax = max(StateL[p],StateR[p])
    h = 1 + sigma*StateL[p]/StateL[r]
    lor = 1/np.sqrt(1 - StateL[v]**2 - StateL[vt]**2)
    A = h*lor*StateL[vt]

    pressure = opt.brentq(lambda pressure: ux(xi,StateL,pressure, sign, A) -
                               Solve_RR_vx(StateL,pressure,sign), pmin, pmax)

    S = StateL[p]/np.power(StateL[r],gamma)
    rho = np.power(pressure/S, 1/gamma)
    h = 1 + sigma*pressure/rho
    Vx = ux(xi,StateL,pressure, sign, A)
    Vt = A*np.sqrt( (1-Vx*Vx)/(h*h +A*A) )
    return rho,Vx,Vt,pressure


if __name__ == "__main__":
    # Test problem
    # StateL = [1.0, 0.0, 0.9, 1000]
    # StateR = [1.0, 0.0, 0.9, .01]

    # StateL = [ 0.999999999999998, 0.0, 0.9, 999.999999999995 ]
    # StateR = [0.999999281059228 ,5.55358611647988e-12, 0.900000124923523, 999.998677473982]

    # print(Solve_RS(StateL, StateR, p))
    # exit(0)

    # StateR = [0.482647065761285, 0.137082381825524, 0.949284868771044,
    #                   296.972419369809]
    # StateL = [0.486091920584473, 0.136066029998986, 0.949050531361542,
    #                   300.513512308036]

    StateR = [0.999999999999559, -4.52552302505939e-17, 0, 999.999999999846]
    StateL = [1.0000009350629, -4.04613917707267e-10, 0, 1000.00074624785]
   # SR Case
    # StateL = [1.0, .5, 0.0, 1]
    # StateR = [.125, 0.0, 0.3, .1]

   # 2R Case
    # StateL = [1.0, 0.0, 0.9, 1]
    # StateR = [.125, 0.5, 0.0, .1]

   # 2S Case
    # StateL = [1.0, 0.5, 0.0, 1]
    # StateR = [.125, 0.0, 0.999, .1]


    Wave3, Wave3Tick = RightWaveType(StateL,StateR)

    # RightShockVT(Wave3Tick,StateR,1)
    # ContactLocation(Wave3)
    # # RightShockVT(StateL,Wave3,-1)
    # print(RareFactionTails(StateL, Wave3,-1))
    # RareState = RarefactionState(.0/.4,StateL,Wave3,-1)

    # print(RareState[vt])
    # print(Wave3[p]/1000)




 # Solve_RR_vx == computeVxb
# class Parent {
# public:
#   double data = 0;
#   void print() { cout << "I am a parent" << endl; }
# };

# class Child : public Parent {
# public:
#   void print() { cout << "I am a child" << endl; }
# };

# int main() {
#   Parent Test;
#   Test.print();
#   Test.data = 10;
#   Child Test2 = *(Child *)&Test;
#   // Child Test2 = *pTest2;
#   Test2.print();
#   cout << Test2.data << endl;

#   return 0;
# }
