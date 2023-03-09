# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 14:05:52 2023

@author: cbri3325
"""


#%% Import functions 

import SimpleITK as sitk
import glob
import os
import pandas as pd
from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from scipy import stats


#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('AIRC24946_R052')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/MRI/baseline/micro/' #Set path to subject directory where cell density map image and GTV image are stored
    out_dir = data_supradir+current+'/MRI/baseline/Prescribed_dose'
              
    cell_map = sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/cellApp_CORRECTED.nii') #Read cellular density map. Units: um-3
    cell_map = cell_map*(10**8) #Convert cell map to Units: 10^4 cm-3
    
    cell_dir = data_supradir+current+'/MRI/baseline/noBone/'
    CTV = sitk.ReadImage(cell_dir+'ctv_voxel_ok_3D.nii.gz')
    CTV = sitk.Cast(CTV, sitk.sitkUInt8)
        
    cell_map = generate_mask(cell_map, CTV)
    cell_map_nda = sitk.GetArrayFromImage(cell_map) #Convert cell map image into a numpy array
    
    #Find indexes of voxels that have values of cellularity
    result = np.where(cell_map_nda>0)
    listOfCoordinates= list(zip(result[0], result[1], result[2])) #Save indexes in a list of coordinates
    
    orig_shape = sitk.GetArrayFromImage(cell_map).shape #Save the shape of the original cell map image
    voxel_volume = (cell_map.GetSpacing()[0]*cell_map.GetSpacing()[1]*cell_map.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
        
    #Define input parameters for TCP modelling: number of cells per voxel (N0) and alpha and beta values per voxel
    CTV_final = CTV
        
    cell_dens = allVoxInt(cell_map, CTV_final) #Define array of values with cellular density>0
    N0 = cell_dens*voxel_volume #Multiply the values of cellular density per voxel by the volume of one voxel to obtain the number of cells per voxel

      
#%%Plot the probability density function vs n cells x vox

    plt.hist(N0, bins=200, density=True)
    plt.xlabel("# Cells per voxel")
    plt.ylabel("Probability density")
    plt.show()
    

#%%# Plot a the kernel density function
    
    import scipy.integrate as integrate
    def integrate_pdf(N):
        return integrate.quad(lambda N: kde(N), 0, N)[0]
    
    def cdf_kde(X):
        "x: list of values of celluarity"
        if np.size(X) > 1:
            return np.array([integrate_pdf(x) for x in X])
        else:
            return integrate_pdf(X)
        
#%%For plotting PDF and CDF
    fig, ax = plt.subplots()
    kde=stats.gaussian_kde(N0)
    X = np.linspace(0, max(N0)*1.1, 200)
    PDF = kde(X)
    CDF = cdf_kde(X)
    ax2 = ax.twinx()
    ax2.fill_between(X, CDF, fc="#53e686", alpha=0.5)
    ax.fill_between(X, PDF, fc="#AAAAFF")
    ax.set_ylabel('Probability density function')
    ax2.set_ylabel('Cumulative density function')
    ax.set_xlabel('N of cells/voxel')
    fig.legend(['PDF', 'CDF'], loc=1)
    plt.show()
    plt.savefig(out_dir+'/PDF_CDF_N_plot.jpg')
    plt.close()
    
    
    def transfer_func(X, Dmin, Dmax, N_dmax):
        """Inputs:
            N cellularity
           Outputs:
            Dose
        """
        return (Dmax-Dmin)*cdf_kde(X)/cdf_kde(N_dmax) + Dmin
    
    Dmin, Dmax = 54, 76
    pdf_vals = tuple(zip(kde(X), X))
    N_dmax = max(pdf_vals)[1]
    dose = transfer_func(N0, Dmin, Dmax, N_dmax)
    
    fig, ax = plt.subplots()
    plt.scatter(N0, dose, marker='x')
    ax2 = ax.twinx()
    ax.plot([N_dmax,N_dmax],[Dmin, max(dose)], 'k--')
    ax.plot([0,max(N0)],[Dmax, Dmax], 'k--')
    ax.annotate('Max Dose {0}'.format(Dmax), (max(N0)*0.7, Dmax*1.02))
    ax.annotate('mode of N density {0}'.format(int(N_dmax)), (N_dmax*1.1, Dmin*1.02))
    ax2.fill_between(X, PDF, fc="#AAAAFF", alpha=0.5)
    ax.set_ylabel('Dose [Gy]')
    ax2.set_ylabel('Probability density function')
    ax.set_xlabel('N of cells/voxel')
    plt.show()
    plt.savefig(out_dir+'/PDF_DOSE_N_plot.jpg')
    plt.close()
    
  #%% Saving the optimised dose values into a prescription dose image
    
    dose_nda = np.ndarray(shape=orig_shape) #Generate a 3D array of the same shape of the cell density map
        
    for d,i in zip(list(dose),listOfCoordinates): #Assign the calculated optimised dose to each voxel, based on the index of the corresponding cellularity voxel
        dose_nda[i]=d
    
    dose_img = sitk.GetImageFromArray(dose_nda) #Convert the numpy array into an image 
    dose_img.CopyInformation(cell_map) #Copy the headers of the cellularity image onto the headers of the optimised dose image
        
    sitk.WriteImage(dose_img, out_dir+'/Dose_optimised_CTV_new.nii')
