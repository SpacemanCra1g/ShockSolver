#!/usr/bin/env python3
import numpy as np
import math
import matplotlib.pyplot as plt
r = 0
v = 1
p = 2
gamma = 1.4


def NewtonsMethod(StateL,StateR,f,fp,guess):
    guess2 = guess - (fside(guess,StateL) + fside(guess,StateR) + (StateR[v] - StateL[v]) )/(fprime(guess,StateL) + fprime(guess,StateR))
    if guess2 < 0.0:
        guess2 = 1.e-10
    
    Cha = abs(guess2 - guess)/( .5*(guess + guess2))

    while Cha > 1e-10:
        guess = guess2
        guess2 = guess - (fside(guess,StateL) + fside(guess,StateR) + (StateR[v] - StateL[v]) )/(fprime(guess,StateL) + fprime(guess,StateR))
        Cha = abs(guess2 - guess)/( .5*(guess + guess2))

    return guess2

def fside(pres,state):
    if (state[p] < pres):
        A = (2/((gamma+1)*state[r]))
        B = ((gamma-1)/(gamma+1))*state[p]
        return (pres - state[p])*math.sqrt(A/(pres +  B) )
    else:
        a = math.sqrt(gamma*state[p]/state[r])
        return (2*a/(gamma-1))*(   (pres/state[p])**((gamma-1)/(2*gamma)) - 1)


def fprime(pres,state):
    if (state[p] < pres):
        A = (2/((gamma+1)*state[r]))
        B = ((gamma-1)/(gamma+1))*state[p]
        return math.sqrt(A/(B + pres))*(1 - (pres-state[p])/(2*(B + pres)))
    else:
        a = math.sqrt(gamma*state[p]/state[r])
        return (1/(state[r]*a))*((pres/state[p])**(-(gamma+1)/(2*gamma)))

def FindPstar(L,R):
    pstar = .5*(L[p] + R[p])
    return NewtonsMethod(L,R,fside,fprime,pstar)

def FindUstar(L,R,pres):
    return ( .5*(L[v] + R[v]) + .5*(fside(pres,R) - fside(pres,L)))

def FindRhoStar(L,R,pres):
    if pres > L[p]:
        pval = pres/L[p]
        gam = (gamma-1)/(gamma+1)
        RhoL = L[r]*( (pval + gam)/(gam*pval + 1))
    else:
        pval = pres/L[p]
        RhoL = L[r]*(pval**(1/gamma))
    
    if pres > R[p]:
        pval = pres/R[p]
        gam = (gamma-1)/(gamma+1)
        RhoR = R[r]*( (pval + gam)/(gam*pval + 1))
    else:
        pval = pres/R[p]
        RhoR = R[r]*(pval**(1/gamma ))

    return RhoL,RhoR

if __name__ == "__main__":
    L = [1,0,1]
    R = [.125, 0, .1]
    pstar = FindPstar(L,R)
    ustar = FindUstar(L,R,pstar)
    rhoLstar, rhoRstar = FindRhoStar(L,R,pstar)

    print(pstar, ustar, rhoLstar, rhoRstar)

    # L = [1, -2, .4]
    # R = [1, 2, .4]
    # pstar = FindPstar(L,R)
    # ustar = FindUstar(L,R,pstar)
    # rhoLstar, rhoRstar = FindRhoStar(L,R,pstar)

    # print(pstar, ustar, rhoLstar, rhoRstar)

    # L = [1,0, 1000]
    # R = [1,0,.01]
    # pstar = FindPstar(L,R)
    # ustar = FindUstar(L,R,pstar)
    # rhoLstar, rhoRstar = FindRhoStar(L,R,pstar)

    # print(pstar, ustar, rhoLstar, rhoRstar)

    # L = [1,0, .01]
    # R = [1,0, 100]
    # pstar = FindPstar(L,R)
    # ustar = FindUstar(L,R,pstar)
    # rhoLstar, rhoRstar = FindRhoStar(L,R,pstar)

    # print(pstar, ustar, rhoLstar, rhoRstar)

    # L = [5.99924,19.5975,460.894]
    # R = [5.99242,-6.19633,46.0950]
    # pstar = FindPstar(L,R)
    # ustar = FindUstar(L,R,pstar)
    # rhoLstar, rhoRstar = FindRhoStar(L,R,pstar)

    # print(pstar, ustar, rhoLstar, rhoRstar)

    time = np.linspace(0,.25,100)
    wave1x = np.linspace(0,.5,50)
    wave1 = abs(wave1x/ustar)
    timestop = [.35]*100

    wave2 = L[v] - math.sqrt(gamma*L[p]/L[r])*math.sqrt( (gamma+1)*pstar/(2*gamma*L[p]) + (gamma-1)/(2*gamma))
    # print(wave2)
    wave2 = abs(wave1x/wave2)

    wave3 = R[v] + math.sqrt(gamma*R[p]/R[r])*math.sqrt( (gamma+1)*pstar/(2*gamma*R[p]) + (gamma-1)/(2*gamma))
    # print(wave2)
    wave3 = abs(wave1x/wave3)

    wave4 = L[v] - math.sqrt(gamma*L[p]/L[r])
    wave4 = abs(wave1x/wave4)

    wave5 = ustar - math.sqrt(gamma*pstar/rhoLstar)
    wave5 = abs(wave1x/wave5)
    wave5 = wave5

    plt.plot(wave1x[:33],wave1[:33])
    plt.plot(-wave1x[:28],wave2[:28],'-', color='black')
    plt.plot(wave1x,wave3)
    plt.plot(-wave1x[:38],wave4[:38],'--', color='black')
    plt.plot(-wave1x[:4],wave5[:4],'--', color='black')
    
    plt.plot(np.linspace(-.5,.5,100),timestop)
    plt.grid()
    plt.show()
    
