# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 02:41:55 2015

@author: JeremyInman
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import csv
from scipy.stats import poisson

def monte_carlo_stock_process(n,Trials):
    print()
    print("Monte Carlo Stock Process and Option Valuation")
    print()
    
    print("Time steps = %s" %n)
    print("Simulations = %s" %Trials)
    print()
    
    #Black-Scholes Pricing Inputs
    S_0=118.79
    K=120.00
    sigma=0.4615
    r=0.0066
    T=240.0/365.0
    h = T/n
    divYield = 0.0000
    
    print("S_0 = %s" %S_0)
    print("K = %s" %K)
    print("sigma = %s" %sigma)
    print("r = %s" %r)
    print("T = %s" %T)
    print("divYield = %s" %divYield)
    print()
    
    
    
    #Allow user to turn jumps on and off
    Jumps_On = True
    
    #used for poisson process
    #remember to create it as a rate (jumps/years)
    #Occurrences of positive asset price jumps, occurrences/years surveyed
    Lambda_up = 5.0/3.0
    Lambda_down = 5.0/3.0
    
    #Mean of jumps surveyed
    alpha_jump_up = 0.1639
    alpha_jump_down = -0.1650
    
    #Standard deviation of jumps surveyed
    sigma_jump_up = 0.0082
    sigma_jump_down = 0.0436
    
    print("Jumps up per year = %s" %Lambda_up)
    print("Jump up mean = %s" %alpha_jump_up)
    print("Jump up sigma = %s" %sigma_jump_up)
    print()
    
    print("Jumps down per year = %s" %Lambda_down)
    print("Jump down mean = %s" %alpha_jump_down)
    print("Jump down sigma = %s" %sigma_jump_down)
    print()
    
    print("Jumps on = %s" %Jumps_On)
    print()
    
    #Poisson probability thresholds
    up_jump1_threshold = 1 - poisson.pmf(1,Lambda_up*h)
    up_jump2_threshold = 1 - poisson.pmf(2,Lambda_up*h)
    up_jump3_threshold = 1 - poisson.pmf(3,Lambda_up*h)
    up_jump4_threshold = 1 - poisson.pmf(4,Lambda_up*h)
    down_jump1_threshold = 1 - poisson.pmf(1,Lambda_down*h)
    down_jump2_threshold = 1 - poisson.pmf(2,Lambda_down*h)
    down_jump3_threshold = 1 - poisson.pmf(3,Lambda_down*h)
    down_jump4_threshold = 1 - poisson.pmf(4,Lambda_down*h)
    
    '''
    print("1 jump up threshold = %s" %up_jump1_threshold)
    print("2 jumps up threshold = %s" %up_jump2_threshold)
    print("3 jumps up threshold = %s" %up_jump3_threshold)
    print("4 jumps up threshold = %s" %up_jump4_threshold)
    print()
    
    print("1 jump down threshold = %s" %down_jump1_threshold)
    print("2 jumps down threshold = %s" %down_jump2_threshold)
    print("3 jumps down threshold = %s" %down_jump3_threshold)
    print("4 jumps down threshold = %s" %down_jump4_threshold)
    print()
    '''
    
    #nname='Monte_Carlo'
    #name=nname+'.csv'

    #file = open(name,'w',newline='')
    #writer = csv.writer(file, quoting=csv.QUOTE_ALL)
    #writer.writerow(['Monte Carlo Simulation']) #file header

    Expiration_Assets = np.zeros(Trials,dtype=np.float32)
    Call_Exercises = np.zeros(Trials,dtype=np.float32)
    Put_Exercises = np.zeros(Trials,dtype=np.float32)
    
    trial=1
    
    Jumps_Up = 0
    Jumps_Down = 0
    
    while trial < Trials + 1:
    
        S=np.zeros(n+1,dtype=np.float32)
        N=np.zeros(n+1,dtype=np.float32)
        D=[]
    
        i=0 #time step counter
    
        for i in range(0,n+1):
            N[i] = i*(T/n)
        
        
        S[0] = S_0
    
        k=1 #time step counter
    
        for k in range(1,n+1):
            
            m_up=0
    
            draw=np.random.random()
            D.append(draw)
    
            if draw > up_jump4_threshold:
                m_up=4
            elif draw > up_jump3_threshold:
                m_up=3
            elif draw > up_jump2_threshold:
                m_up=2
            elif draw > up_jump1_threshold:
                m_up=1
            else:
                m_up=0
            
            Jumps_Up += m_up #keeping a running total of jumps up in all trials
            
            m_down=0
    
            draw=np.random.random()
    
            if draw > down_jump4_threshold:
                m_down=4
            elif draw > down_jump3_threshold:
                m_down=3
            elif draw > down_jump2_threshold:
                m_down=2
            elif draw > down_jump1_threshold:
                m_down=1
            else:
                m_down=0
            
            Jumps_Down += m_down #keeping a running total of jumps down in all trials
            
            
            if Jumps_On == True:
                
                jumpAdjustment_up = Lambda_up*(math.exp(alpha_jump_up) - 1)
                jumpAdjustment_down = Lambda_down*(math.exp(alpha_jump_down) - 1)
                
                Jump_Drift_Up = math.exp(m_up*(alpha_jump_up - 0.5*sigma_jump_up**2))
                Jump_Drift_Down = math.exp(m_down*(alpha_jump_down - 0.5*sigma_jump_down**2))
                
                Jump_Noise_Rand_Sum = 0
                j=1
                
                while j < m_up+1:
                    Jump_Noise_Rand_Sum += np.random.normal(0,1)
                    j += 1
                
                Jump_Noise_Up = math.exp(sigma_jump_up*Jump_Noise_Rand_Sum)
                
                Jump_Noise_Rand_Sum = 0
                j=1
                
                while j < m_down+1:
                    Jump_Noise_Rand_Sum += np.random.normal(0,1)
                    j += 1
                
                Jump_Noise_Down = math.exp(sigma_jump_down*Jump_Noise_Rand_Sum)
                
            else:
                jumpAdjustment_up = 0
                jumpAdjustment_down = 0
                Jump_Drift_Up = 1
                Jump_Noise_Up = 1
                Jump_Drift_Down = 1
                Jump_Noise_Down = 1
            
            Drift = math.exp((r - jumpAdjustment_up - jumpAdjustment_down - divYield - 0.5*sigma**2)*h)
            Noise = math.exp(sigma*math.sqrt(h)*np.random.normal(0,1))
            
            S[k] = S[k-1]*Drift*Noise*Jump_Drift_Up*Jump_Noise_Up*Jump_Drift_Down*Jump_Noise_Down
            
            if k == n:
                Expiration_Assets[trial-1] = S[k]

        #file = open(name,'a',newline='') #open the file to append
        #writer = csv.writer(file, quoting=csv.QUOTE_ALL)
        #writer.writerow(['Trial', trial]) #adjusts because the counter starts at 0
    
        #for values in N,S: #write the trial's values
                    #writer.writerow(values)
        
        plt.plot(N,S)
        trial+=1
    
    #verify mean of uniform draws, should be 0.5
    mean=np.mean(D)
    print("uniform draw mean= %s" %mean)
    print("Jumps Up = %s" %Jumps_Up)
    print("Jumps Down = %s" %Jumps_Down)
    print()
    
    Average_Expiration_Asset_Price = np.mean(Expiration_Assets)
    
    print("Average Asset Price = %s" %Average_Expiration_Asset_Price)
    print()
    
    i=0
    
    for i in range(0,Trials):
            Call_Exercises[i] = max(Expiration_Assets[i] - K,0)
            Put_Exercises[i] = max(K - Expiration_Assets[i],0)
    
    
    Call_Price = math.exp(-r*T)*np.mean(Call_Exercises)
    Put_Price = math.exp(-r*T)*np.mean(Put_Exercises)
    
    Call_Price_STDev = np.std(Call_Exercises)*(1/math.sqrt(Trials))
    Put_Price_STDev = np.std(Put_Exercises)*(1/math.sqrt(Trials))
    
    print("Call Price = %s" %Call_Price)
    print("Call Price STDev = %s" %Call_Price_STDev)
    print()
    
    print("Put Price = %s" %Put_Price)
    print("Put Price STDev = %s" %Put_Price_STDev)
    print()
    
# main program starts here
monte_carlo_stock_process(500,2000)