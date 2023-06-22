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

def lc(x,a,b):
    
    '''Find LC given D_98%
    
    x = dose delivered to 98% of GTV for each patient (array)'''
    return 1/(1+np.exp(a+b*x))

#%%Import DVH info dataframe
dvh_df = pd.read_excel('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/DVH_D95_RayStation.xlsx', sheet_name='Sheet1')

D95 = dvh_df['D95'].to_numpy()
LC_true = dvh_df['LC'].to_numpy()

#%%D95
#%% Fit the data
popt,pcov=curve_fit(lc, D95, LC_true, p0=[1000, -20])
perr = np.sqrt(np.diag(pcov))

#%% Evaluate goodness of the fit

LC_pred = lc(D95,*popt)
Fit_QA= r2_score(LC_true, LC_pred)

#%%Plot results

# xmodel = np.arange(math.floor(np.max(D98))+2)
xmodel = np.arange(100)
ymodel = lc(xmodel, *popt)

plt.plot(D95, LC_true, "go", label="True LC")
plt.plot(xmodel, ymodel, "g--", label="Fit - R^2="+str(round(Fit_QA, 2)))
plt.xlabel('D95% Dose [Gy]')
plt.ylabel('LC')
plt.axes().xaxis.set_minor_locator(MultipleLocator(0.5))
plt.legend()
plt.savefig('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/LC_fit_D95_RayStation_sigmoid.jpg')
plt.show()

print('For LC vs D95 a was ',popt[0],' and b was ',popt[1])
