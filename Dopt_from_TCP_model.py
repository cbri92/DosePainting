# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 08:13:35 2022

@author: cbri3325
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def theta(x, N0, alpha):
    """Optimisation problem for finding db array that minimizes the
    problem where TCP = 1, for known N0 and alpha arrays
    """
    prod = np.prod( np.exp( - N0 * np.exp(-alpha*x) )) #This model uses a Poisson distribution with the linear quadratic model to predict the probability of tumor control
    penalty = x*0.1/100
    # penalty = 0
    return 1 - prod + penalty

# define number of voxels and randomly assign alpha and N0
N = 3
alpha = np.random.normal(0.30, 0.05, size=(N,N,N)).flatten()
N0 = np.random.normal(loc=7400, scale=1000, size=(N,N,N)).flatten()

# for graphing below only
bnds =(0,200) #Range of values of x that the minimization function will search from
x_vals = np.arange(0,200)
vals=[]
for x in np.arange(0,200):
    vals.append(np.prod( np.exp( - N0 * np.exp(-alpha*x) )))

# optimisation problem
x0 = 20
bnds = [bnds]
x = minimize(theta, x0, args=(N0,alpha), method='Powell', tol=1e-10, bounds = bnds)

val = np.array(vals)
plt.plot(x_vals,val)
plt.plot([x.x,x.x],[0,1])
plt.xlabel("Dose (Gy)")
plt.ylabel("TCP")

plt.plot(x_vals[:-1],np.diff(val))
plt.plot([x.x,x.x],[0,max(np.diff(val))])
plt.xlabel("Dose (Gy)")
plt.ylabel("TCP step difference")

# way of interpolating dose plan
# metric of severity of voxel: increasing cell number N0 and decreasing alpha is higher risk
param = (N0/alpha)
param / np.sum(param), param # weighting of each voxel

# Normalise [0,1] and assign 50 Gy to 0 (least at risk) to optimisation number previously determined
# User linear interpolation to assign all other values for heterogenous plan
norm01 = (param - np.min(param))/np.ptp(param)
norm01

def dose_calc(arr):
    LL, UL = 50,x.x
    return (UL - LL)*arr + LL

# heterogenous dose plan
dose_calc(norm01)

# Compare new probabilty
np.prod( np.exp( - N0 * np.exp(-alpha*dose_calc(norm01)) ))

# To old probability with optimised dose to entire brain
np.prod( np.exp( - N0 * np.exp(-alpha*x.x)) )
