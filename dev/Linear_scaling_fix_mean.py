# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:49:29 2023

@author: cbri3325
"""

import numpy as np
import matplotlib.pyplot as plt

a = np.random.randint(50, size=20)
a_max = max(a)
a_min = min(a)
a_mean = np.mean(a)


plt.hist(a, bins=50)
plt.xlabel("Value")
plt.ylabel("Counts")
plt.show()


d_max = (((a_max-a_min)*(74-68))/(a_mean-a_min))+68

m = ((74-68)/(a_mean-a_min))

y = m*a + 68

plt.hist(y, bins=50)
plt.xlabel("dose")
plt.ylabel("Counts")
plt.show()

y_mean = np.mean(y)
