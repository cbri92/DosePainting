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
from sklearn.metrics import explained_variance_score


#%%Define LC function

# def lc(x, N0):
    
#     '''Find LC given D_98%
    
#     alphaPOP = 0.36 Gy-1 from fit of data (float)
#     N0 = Number of cells in GTV - parameter to be fund (array)
#     x = dose delivered to 98% of GTV for each patient (array)'''
#     return np.exp( - N0 * np.exp(-alphaPOP*x))

# alphaPOP=0.3626439540025473 #Gy-1

def lc(x, a):
    
    '''Find LC given D_98%
    
    a = alpha Gy-1 to be found (float)
    10**7 = Number of cells in GTV - from literature
    x = dose delivered to 98% of GTV for each patient (array)'''
    return np.exp( - (10**7) * np.exp(-a*x))

#%%Import DVH info dataframe
dvh_df = pd.read_excel('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/DVH_D95_RayStation.xlsx', sheet_name='Sheet1')

# N0 = dvh_df['N0'].to_numpy()
# N0=10**7

D95 = dvh_df['D95'].to_numpy()
LC_true = dvh_df['LC'].to_numpy()

#%%D95
#%% Fit the data

a, pcov = curve_fit(lc, D95, LC_true, p0=[0.36], bounds=(0, np.inf))
# N0, pcov = curve_fit(lc, D95, LC_true, p0=[10**7], bounds=(0, np.inf))
perr = np.sqrt(np.diag(pcov))

#%% Evaluate goodness of the fit

LC_pred = lc(D95, a)
# LC_pred = lc(D95, N0)
Fit_QA= r2_score(LC_true, LC_pred)
# Fit_QA2= explained_variance_score(LC_true, LC_pred)

#%%Plot results

# xmodel = np.arange(math.floor(np.max(D98))+2)
xmodel = np.arange(100)
# ymodel = lc(xmodel, N0)
ymodel = lc(xmodel, a)

plt.plot(D95, LC_true, "go", label="True LC")
plt.plot(xmodel, ymodel, "g--", label="Fit - R^2="+str(round(Fit_QA, 2)))
plt.xlabel('D95% Dose [Gy]')
plt.ylabel('LC')
plt.axes().xaxis.set_minor_locator(MultipleLocator(0.5))
plt.legend()
plt.savefig('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/LC_fit_D95_RayStation.jpg')
plt.show()

print('Alpha for for LC vs D95 from RayStation was ',a)