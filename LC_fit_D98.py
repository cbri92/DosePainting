# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:03:06 2023

@author:  Caterina Brighi
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
import math
from sklearn.metrics import r2_score


#%%Define LC function

def lc(x, N0):
    
    '''Find LC given D_98%
    
    alphaPOP = 0.19 Gy-1 from fit of data (float)
    N0 = Number of cells in GTV - parameter to be fund (array)
    x = dose delivered to 98% of GTV for each patient (array)'''
    return np.exp( - N0 * np.exp(-alphaPOP*x))

alphaPOP=0.3626439540025473 #Gy-1

#%%Import DVH info dataframe
dvh_df = pd.read_excel('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/DVH_info.xlsx', sheet_name='DVH stats')

# N0 = dvh_df['N0'].to_numpy()
# N0=10**7
D98 = dvh_df['D_98% [Gy]'].to_numpy()
LC_true = dvh_df['LC'].to_numpy()

#%% Fit the data

N0, pcov = curve_fit(lc, D98, LC_true, p0=[0.05], bounds=(0, np.inf))
perr = np.sqrt(np.diag(pcov))

#%% Evaluate goodness of the fit

LC_pred = lc(D98, N0)
Fit_QA= r2_score(LC_true, LC_pred)

#%%Plot results

# xmodel = np.arange(math.floor(np.max(D98))+2)
xmodel = np.arange(100)
ymodel = lc(xmodel, N0)

plt.plot(D98, LC_true, "ro", label="True LC")
plt.plot(xmodel, ymodel, "r--", label="Fit - R^2="+str(round(Fit_QA, 2)))
plt.xlabel('D98 Dose [Gy]')
plt.ylabel('LC')
plt.axes().xaxis.set_minor_locator(MultipleLocator(0.5))
plt.legend()
plt.savefig('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/LC_fit_D98.jpg')
plt.show()


