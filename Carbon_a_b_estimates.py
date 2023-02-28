# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 08:27:06 2023

@author: Caterina Brighi
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
import math
from statsmodels.stats.weightstats import DescrStatsW

def lq(x, a, b):
    return -(a*x+b*x*x)


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))

x_30 = np.array([0,1,2,3,4])
y_30 = np.array([1,0.39,0.195,0.073,0.0185])

x_50 = np.array([0,1,2,3,4])
y_50 = np.array([1,0.285,0.088,0.0415,0.0165])

x_70 = np.array([0,0.5,1,2,3])
y_70 = np.array([1,0.735,0.34,0.075,0.019])

result_dir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/'
Results = pd.ExcelWriter(result_dir +'C_curveFitting_results.xlsx')

#%%Method 1

res_df = pd.DataFrame(columns=['Radiation type', 'alpha', 'alpha std', 'beta', 'beta std', 'Weighted alpha', 'Weighted alpha std', 'Weighted beta', 'Weighted beta std'])

(a30,b30), pcov30 = curve_fit(lq, x_30, np.log(y_30), p0=[0.05,0.05], bounds=(0, np.inf))
perr30 = np.sqrt(np.diag(pcov30))
res30 = {'Radiation type':'Carbon 30', 'alpha':a30, 'alpha std':perr30[0], 'beta':b30, 'beta std':perr30[1]}

(a50,b50), pcov50 = curve_fit(lq, x_50, np.log(y_50), p0=[0.05,0.05], bounds=(0, np.inf))
perr50 = np.sqrt(np.diag(pcov50))
res50 = {'Radiation type':'Carbon 50', 'alpha':a50, 'alpha std':perr50[0], 'beta':b50, 'beta std':perr50[1]}

(a70,b70), pcov70 = curve_fit(lq, x_70, np.log(y_70), p0=[0.05,0.05], bounds=(0, np.inf))
perr70 = np.sqrt(np.diag(pcov70))
res70 = {'Radiation type':'Carbon 70', 'alpha':a70, 'alpha std':perr70[0], 'beta':b70, 'beta std':perr70[1]}

#Append results to dataframe
res_df = res_df.append(res30, ignore_index=True)
res_df = res_df.append(res50, ignore_index=True)
res_df = res_df.append(res70, ignore_index=True)

#Calculate weighted mean of alpha and beta
weig_alpha, std_wAlpha = weighted_avg_and_std(res_df['alpha'], res_df['alpha std']) 
weig_beta, std_wBeta = weighted_avg_and_std(res_df['beta'], res_df['beta std'])
weighted_values ={'Weighted alpha':weig_alpha, 'Weighted alpha std':std_wAlpha, 'Weighted beta':weig_beta, 'Weighted beta std':std_wBeta}

res_df = res_df.append(weighted_values, ignore_index=True)

res_df.to_excel(Results, sheet_name='CurveFit', index=False)

#Plot results of the fit
xmodel = np.arange(6)
ymodel_30 = np.exp(lq(xmodel, a30, b30))  # Note: take exp() as inverse log()
ymodel_50 = np.exp(lq(xmodel, a50, b50))
ymodel_70 = np.exp(lq(xmodel, a70, b70))

plt.plot(x_30, y_30, "ro", label="Carbon 30 Experiment")
plt.plot(xmodel, ymodel_30, "r--", label="Carbon 30 Model")
plt.plot(x_50, y_50, "bo", label="Carbon 50 Experiment")
plt.plot(xmodel, ymodel_50, "b--", label="Carbon 50 Model")
plt.plot(x_70, y_70, "go", label="Carbon 70 Experiment")
plt.plot(xmodel, ymodel_70, "g--", label="Carbon 70 Model")
plt.yscale('log')
plt.xlim(0, 6)
plt.ylim(0.01,1)
plt.xlabel('Dose [Gy]')
plt.ylabel('Survival Fraction')
plt.axes().xaxis.set_minor_locator(MultipleLocator(0.5))
plt.legend()
plt.savefig(result_dir+'Carbon_curveFit.jpg')
plt.show()

#%%Method 2

from lmfit import Model
import pybroom

model = Model(lq)
model.set_param_hint('a', min=0)
model.set_param_hint('b', min=0)
params = model.make_params(a=0.05, b=0.05)

Res_df = pd.DataFrame(columns=['Radiation type', 'alpha', 'alpha std', 'beta', 'beta std', 'Weighted alpha', 'Weighted alpha std', 'Weighted beta', 'Weighted beta std'])

result30 = model.fit(np.log(y_30), params, x=x_30)
print('Results Carbon 30')
print(result30.fit_report())
df30=pybroom.tidy_lmfit_result(result30)
df30=df30.set_index('name')
Res30 = {'Radiation type':'Carbon 30', 'alpha':df30._get_value('a', 'value'), 'alpha std':df30._get_value('a', 'stderr'), 'beta':df30._get_value('b', 'value'), 'beta std':df30._get_value('b', 'stderr')}

result50 = model.fit(np.log(y_50), params, x=x_50)
print('Results Carbon 50')
print(result50.fit_report())
df50=pybroom.tidy_lmfit_result(result50)
df50=df50.set_index('name')
Res50 = {'Radiation type':'Carbon 50', 'alpha':df50._get_value('a', 'value'), 'alpha std':df50._get_value('a', 'stderr'), 'beta':df50._get_value('b', 'value'), 'beta std':df50._get_value('b', 'stderr')}

result70 = model.fit(np.log(y_70), params, x=x_70)
print('Results Carbon 70')
print(result70.fit_report())
df70=pybroom.tidy_lmfit_result(result70)
df70=df70.set_index('name')
Res70 = {'Radiation type':'Carbon 70', 'alpha':df70._get_value('a', 'value'), 'alpha std':df70._get_value('a', 'stderr'), 'beta':df70._get_value('b', 'value'), 'beta std':df70._get_value('b', 'stderr')}

#Append results to dataframe
Res_df = Res_df.append(Res30, ignore_index=True)
Res_df = Res_df.append(Res50, ignore_index=True)
Res_df = Res_df.append(Res70, ignore_index=True)

#Calculate weighted mean of alpha and beta
Weig_alpha, Std_wAlpha = weighted_avg_and_std(Res_df['alpha'], Res_df['alpha std']) 
Weig_beta, Std_wBeta = weighted_avg_and_std(Res_df['beta'], Res_df['beta std'])
Weighted_values ={'Weighted alpha':Weig_alpha, 'Weighted alpha std':Std_wAlpha, 'Weighted beta':Weig_beta, 'Weighted beta std':Std_wBeta}

Res_df = Res_df.append(Weighted_values, ignore_index=True)

Res_df.to_excel(Results, sheet_name='lmfit', index=False)

#Plot results of the fit
xmodel = np.arange(6)
ymodel_30 = np.exp(result30.eval(x=xmodel))
ymodel_50 = np.exp(result50.eval(x=xmodel))
ymodel_70 = np.exp(result70.eval(x=xmodel))

plt.plot(x_30, y_30, "ro", label="Carbon 30 Experiment")
plt.plot(xmodel, ymodel_30, "r--", label="Carbon 30 Model")
plt.plot(x_50, y_50, "bo", label="Carbon 50 Experiment")
plt.plot(xmodel, ymodel_50, "b--", label="Carbon 50 Model")
plt.plot(x_70, y_70, "go", label="Carbon 70 Experiment")
plt.plot(xmodel, ymodel_70, "g--", label="Carbon 70 Model")
plt.yscale('log')
plt.xlim(0, 6)
plt.ylim(0.01,1)
plt.xlabel('Dose [Gy]')
plt.ylabel('Survival Fraction')
plt.axes().xaxis.set_minor_locator(MultipleLocator(0.5))
plt.legend()
plt.savefig(result_dir+'Carbon_lmfit.jpg')
plt.show()

#%%Save results

Results.save()