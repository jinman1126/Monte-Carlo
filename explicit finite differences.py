# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 02:41:55 2015

@author: JeremyInman
"""


import numpy as np
import matplotlib.pyplot as plt
import math

def bs_finite_exp():
    print()
    print('Black-Scholes PDE, Numerical Estimation, Explicit Finite Differences')
    print()
#Set the current stock price S curr, strike price K, volatility sigma, risk-free rate r, time toexpiration T
    S_curr = 105.99
    K = 105.00
    sigma = 0.3128
    r = 0.0017
    T = 351.0/365.0
    divYield = 0.017
#Set number of asset price steps J
    J = 100
#Set number of time steps M and time increment deltaT
    M = 1000
#Set size of deltaS and deltaT
    S_max = 2*S_curr
    deltaS = S_max/J
    deltaT = T/M
        
    print('Strike Price K = %s' %K)
    print('Volatility sigma = %s' %sigma)
    print('Risk-Free Rate r = %s' %r)
    print('Time to expiration T = %s' %T)
    print('Dividend yield divYeild = %s' %divYield)
    print()
    
    print('Asset price steps J = %s' %J)
    print('Time increment steps M = %s' %M)
    print('Size of Asset price steps = %s' %deltaS)
    print('Size of the time steps = %s' %deltaT)
    
        
#Create (J-1) X (J-1) Matrix A
#Create column vector of asset prices S
#Create column vector of option prices C hat
    A=np.zeros([J-1,J-1],dtype=float)
    S=np.zeros(J-1,dtype=float)
    C_hat=np.zeros(J-1,dtype=float)
        
    S[0] = deltaS
    j=0 #row counter
    for j in range(1,J-1):
        S[j] = S[j-1] + deltaS
        
    j=0 #row counter
    for j in range(0,J-1):
        C_hat[j] = max(S[j]-K,0)
            
    j=0 #row counter
    k=0 #column counter
    for j in range (0,J-1):
        for k in range (0,J-1):
            if k == j:
                A[j,j] = 1 - (sigma**2)*((j+1)**2)*deltaT - r*deltaT
            elif k == j-1:
                A[j,k] = 0.5 * ((sigma**2)*((j+1)**2)*deltaT - (r - divYield)*(j+1)*deltaT)
            elif k == j+1:
                A[j,k] = 0.5 * ((sigma**2)*((j+1)**2)*deltaT + (r - divYield)*(j+1)*deltaT)
                    
    C_0 = 0
    C_start = C_hat
    m=1 #time increment counter
    while m <= M:
        #Boundary values calculated separately
        #Must incorporate lowest and highest option prices, which are not in matrix
        C_min = 0.5*((sigma**2)*((1)**2)*deltaT - (r - divYield)*(1)*deltaT)*C_0 + (1 - (sigma**2)*((1)**2)*deltaT - r*deltaT)*C_hat[0] + 0.5*((sigma**2)*((1)**2)*deltaT + (r - divYield)*(1)*deltaT)*C_hat[1]
        C_max = 0.5*((sigma**2)*((J-1)**2)*deltaT - (r - divYield)*(J-1)*deltaT)*C_hat[J-3] + (1 - (sigma**2)*((J-1)**2)*deltaT - r*deltaT)*C_hat[J-2] + 0.5*((sigma**2)*((J-1)**2)*deltaT + (r - divYield)*(J-1)*deltaT)*(S_max - K*math.exp(- r*m*deltaT))
        #Perform matrix multiplication Cm+1 = ACm
        C_hat = A.dot(C_hat)
        #Update boundary values
        C_hat[0] = C_min
        C_hat[J-2] = C_max
        #Capture option prices during backward walk for graphing
        if m == M/4:
            C1 = C_hat
        elif m == 2*M/4:
            C2 = C_hat
        elif m == 3*M/4:
            C3 = C_hat
        elif m == 4*M/4:
            C4 = C_hat
        m=m+1
    J_2 = float((J-2)/2)
    option_price = float(C_hat[J_2])    
    print('Option Price = %s' %option_price)
    plt.plot(S,C_start,'b--',S,C1,'r--',S,C2,'g--',S,C3,'r--',S,C4,'b--')
    # main program starts here
bs_finite_exp()