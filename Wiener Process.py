"""
Created on Wed Sep 28 17:38:13 2016

@author: jeremy
"""

#Wiener Process
#Spec_Readings

import random
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import csv

def wiener_process():
    
    #lets get some of the inputs we need
    print('I recommend starting with 10 trials at 200-500 simulations, adjust from there.')
    print()
    print('I have not put exception handling in here yet, so only positive integers for the trials and simulations!')
    print()
    print('As for the file name, exclude any special characters and spaces')
    print()
    t= int(input(print('How many trials would you like to run?')))
    n = int(input(print('How many simulations will we be running in each trial?')))
    
    #now we have to create our output file
    nname=input('please name your file (without any extension tags):  ') #i add the .csv below
    name=nname+'.csv'  
    file = open(name,'w',newline='')
    writer = csv.writer(file, quoting=csv.QUOTE_ALL)
    writer.writerow(['Wiener Process']) #just a header for the file
    
    print()
    print('Wiener Process')
    print()
    print('Calculating...')
    
    
    #t=2 #trial number
    tc=1 #trial counter
    
    while tc<t+1:
        
        
        #n=1000
        print('trial number: ', tc) #give them something to look at while it is loading
        #print('n=%s'%n) #return value as confirmation
    
        #create our data sets
        S=np.zeros(n+1,dtype=float) #will give result of last S calculation
        N=np.zeros(n+1,dtype=float)
        
        for i in range(0,n+1):
        
            N[i]=i
    
        #print('N=%s'%N)

        k=1 #row counter
            
        for k in range(1,n+1):
                
            i=0
            sum_i=0
            
            if k==n: #here, we are going to capture the S matrix when it is complete
                
                for i in range(0,k):
                        
                    sum_i=sum_i + np.random.normal(0,1,1)
                        
                S[k]=(1/math.sqrt(n))*sum_i #calculate the value
                tc+=1 #add to the trial counter
                plt.plot(N,S) #plot the trial
                file = open(name,'a',newline='') #open the file to append
                writer = csv.writer(file, quoting=csv.QUOTE_ALL)
                writer.writerow(['Trial', tc-1]) #adjusts because the counter starts at 0
    
    
                for values in N,S: #write the trial's values
                    writer.writerow(values)
        
                file.close()
                
            else: #for most cases, we will just run this until we get to k=n
                for i in range(0,k):
                        
                    sum_i=sum_i + np.random.normal(0,1,1)
                    
                S[k]=(1/math.sqrt(n))*sum_i
                   
    print('S=%s'%S)
     
    #starts here
    
wiener_process()
