# -*- coding: utf-8 -*-
"""
Created on Tue May  2 18:19:20 2023

@author: Caterina Brighi

This script finds alphaPOP through an optimization problem, 
which consisted in finding the value of α (αPOP) for Np patients so that the population TCP (TCPPOP), 
modelled through a Poissonian statistics, matched the observed LCPOP. 

Following method published in supplementary material Step 3 of: Buizza et al. Radiotherapy and Oncology 137 (2019) 32–37.
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt

#%%Define constants

TCPpop = 0.48 #LC of SC patients 24/50 with LC

#%%Define function for optimization problem

def function(alphaPOP, Np, TCPpop, N0, Dlq):
    
    '''Optimization problem for finding alphaPOP that minimizes the
    problem where TCP = 0.48, for known N0 and Dlq arrays.
    
    alphaPOP = variable to be found (float)
    Np = number of patients (int)
    TCPpop = LC of population (float)
    N0 = Number of cells in GTV for each patient (array)
    Dlq = dose lq for each patient (array)'''
    
    somma = np.sum( np.exp( - N0 * np.exp(-alphaPOP*Dlq) ))
    
    # somma = 0
    # for i in range(len(N0)):
    #     print(str(i))
    #     somma += np.exp( - N0[i] * np.exp(-alphaPOP*Dlq[i]) )
        
    f = -Np*TCPpop+somma
    return f

#%%Import DVH info dataframe
dvh_df = pd.read_excel('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/DVH_info.xlsx', sheet_name='DVH stats')

N0 = dvh_df['N0'].to_numpy()
Dlq = dvh_df['Dlq [Gy]'].to_numpy()


#%% Find alphaPOP through optimization function

Np = len(N0)

alphaPOP0 = 0.2
x = minimize(function, alphaPOP0, args=(Np, TCPpop, N0, Dlq), method='Nelder-Mead', tol=1e-20)
alphaPOP=x.x[0]


alphaPOP0_vals = np.linspace(0,2,100)
opt = [function(alpha, Np, TCPpop, N0, Dlq) for alpha in alphaPOP0_vals]

plt.plot(alphaPOP0_vals, opt)
plt.plot([x.x,x.x],[min(opt),max(opt)])
plt.legend(['f','alphaPOP: '+str(round(x.x[0],2))], loc ="lower right")
plt.xlabel("alpha (Gy-1)")
plt.ylabel("f")
plt.show()
