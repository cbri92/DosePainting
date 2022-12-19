# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:31:51 2022

@author: Caterina Brighi
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


def theta(x, N0, alpha, beta, n):
    """Optimisation problem for finding db array that minimizes the
    problem where TCP = 1, for the following known arrays:
        N0: cellularity
        alpha: radiosensitivity parameter
        beta: radiosensitivity parameter
        n: number of fractions
    """
    prod = np.prod( np.exp( - N0 * np.exp(-alpha*x-(beta*(x**2))/n) )) #This model uses a Poisson distribution with the linear quadratic model to predict the probability of tumor control
    return 1 - prod


def dose_calc(array, LL, UL):
    """This function calculates the optimised dose given:
        - an array with the per voxel metric of severity calculated from the minimization problem;
        - LL : the lower dose value that should be delivered;
        - UL : the maximum dose value that cannot be exceeded.
    """
    return (UL - LL)*array + LL


#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/MRI/baseline/micro/' #Set path to subject directory where cell density map image and GTV image are stored
    if not os.path.exists(subj_dir+'Optimised_dose'):#if it does not already exist, create a directory where the optimized dose will be saved
        os.mkdir(subj_dir+'Optimised_dose')
    out_dir = subj_dir+'Optimised_dose'
    
    #Unzip cellApp_CORRECTED.nii.gz file and delete original gzipped file
    for filename in glob.glob(subj_dir +'cellApp_CORRECTED.nii.gz'):
        gunzip_shutil(filename, filename[:-3])
 
    print('Optimizing the dose prescription for '+current) 
          
    cell_map = sitk.ReadImage(subj_dir+'cellApp_CORRECTED.nii') #Read cellular density map. Units: um-3
    cell_map = cell_map*(10**8) #Convert cell map to Units: 10^4 cm-3
    # cell_map = cell_map*(10**12) #Convert cell map to Units: cm-3
    
    alpha_map = sitk.ReadImage(subj_dir+'alpha_inDWI.nii') 
    beta_map = sitk.ReadImage(subj_dir+'beta_inDWI.nii') 
    
    CTV = sitk.ReadImage(subj_dir+'CTVonDWI.nii') #Read CTV
    # CTV.SetOrigin(cell_map.GetOrigin())
    
    cell_map_nda = sitk.GetArrayFromImage(cell_map) #Convert cell map image into a numpy array
    # alpha_map_nda = sitk.GetArrayFromImage(alpha_map) #Convert alpha map image into a numpy array
    # beta_map_nda = sitk.GetArrayFromImage(beta_map) #Convert beta map image into a numpy array    
    
    #Find indexes of voxels that have values of alpha
    result = np.where(cell_map_nda>0)
    listOfCoordinates= list(zip(result[0], result[1], result[2])) #Save indexes in a list of coordinates
    # for cord in listOfCoordinates:
    #     print(cord)    

    orig_shape = sitk.GetArrayFromImage(cell_map).shape #Save the shape of the original cell map image
    voxel_volume = (cell_map.GetSpacing()[0]*cell_map.GetSpacing()[1]*cell_map.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
    
    #Define input parameters for TCP modelling: number of cells per voxel (N0) and alpha and beta values per voxel
    CTV_final = (cell_map>0) & (alpha_map>0) #Create a mask where both alpha and cellularity values are >0
    
    cell_dens = allVoxInt(cell_map, CTV_final) #Define array of values with cellular density>0
    N0 = cell_dens*voxel_volume #Multiply the values of cellular density per voxel by the volume of one voxel to obtain the number of cells per voxel
    
    alpha = allVoxInt(alpha_map, CTV_final) #Define array of values with alpha
    beta = allVoxInt(beta_map, CTV_final) #Define array of values with beta
    
    # alpha = np.random.normal(0.30, 0.05, size=len(N0)) #Define array of values of alpha (radiosensitivity parameter)
    # beta = np.random.normal(0.15, 0.01, size=len(N0)) #Define array of values of beta (radiosensitivity parameter)
    
    n=16 #Number of fractions

#%%Optimization
    
    # for graphing below only
    x_vals = np.arange(0,200)
    #TCP Poisson model
    vals=[np.prod( np.exp( - N0 * np.exp(-alpha*x-(beta*(x**2))/n) )) for x in x_vals]
    
    TCP = np.array(vals)
    VAL_max = np.isclose(TCP, 0.99, rtol=1e-02)
    MAX_DOSE = x_vals[VAL_max][0]
    print('MAX DOSE', MAX_DOSE)
    MIN_DOSE = x_vals[TCP>0.01][0]
    print('MIN DOSE', MIN_DOSE)
    
    #Plot the TCP vs Dose
    plt.plot(x_vals,TCP)
    plt.plot([MAX_DOSE,MAX_DOSE],[0,1])
    plt.plot([MIN_DOSE,MIN_DOSE],[0,1])
    plt.legend(['TCP','Max dose value: '+str(MAX_DOSE), 'Min dose value: '+str(MIN_DOSE)], loc ="upper left")
    plt.xlabel("Dose (Gy)")
    plt.ylabel("TCP")
    plt.savefig(out_dir+'/TCP_dose_plot.jpg')
    plt.show()

    #Plot the probability density function vs dose
    plt.plot(x_vals[:-1],np.diff(TCP))
    plt.plot([MAX_DOSE,MAX_DOSE],[0,max(np.diff(TCP))])
    plt.plot([MIN_DOSE,MIN_DOSE],[0,max(np.diff(TCP))])
    plt.legend(['PDF','Max dose value: '+str(MAX_DOSE), 'Min dose value: '+str(MIN_DOSE)], loc ="upper left")
    plt.xlabel("Dose (Gy)")
    plt.ylabel("Probability density")
    plt.savefig(out_dir+'/PDF_dose_plot.jpg')
    plt.show()
    
#%% Way of interpolating dose to obtain heterogeneous dose plan, such that to target higher dose to more radioresistant regions, whilst sparing dose to radiosensitive tissue (and surrounding healthy tissues)

    # Define a metric of severity of voxel: increasing cell number N0 and decreasing alpha is higher risk
    param = (N0/(alpha)) #Array containing metric of severity estimated per voxel
    param / np.sum(param), param # weighting of each voxel
    
    # Normalise [0,1] and assign 50 Gy (minimum dose to be delivered to the voxel in the GTV) to 0 (least at risk) to optimisation number previously determined
    norm01 = (param - np.min(param))/np.ptp(param)
    norm01
    
    # Develop heterogenous dose plan by using linear interpolation to assign all other values
    Dopt = dose_calc(norm01,MIN_DOSE,MAX_DOSE)
    
    # Compare new TCP with TCP that would be obtained by deliverying an homogeneous dose boost to the GTV to maximum dose estimated in the optimization
    TCP_heter =np.prod( np.exp( - N0 * np.exp(-alpha*Dopt-(beta*Dopt**2)/n))) #TCP obtained with optimised heterogeneous dose distribution
    TCP_boost = np.prod( np.exp( - N0 * np.exp(-alpha*MAX_DOSE-(beta*MAX_DOSE**2)/n))) # TCP with optimised homogeneous max dose boost to entire GTV
    print('TCP optimised heterogeneous dose plan: ',TCP_heter)
    print('TCP optimised homogeneous max dose boost plan: ',TCP_boost)
    
    plt.hist(Dopt, bins=20)
    plt.xlabel("Dose (Gy)")
    plt.ylabel("Counts")
    plt.savefig(out_dir+'/Dose_histo.jpg')
    plt.close()
    
#%% Saving the optimised dose values into a prescription dose image

    dose_nda = np.ndarray(shape=orig_shape) #Generate a 3D array of the same shape of the cell density map
    
    for d,i in zip(list(Dopt),listOfCoordinates): #Assign the calculated optimised dose to each voxel, based on the index of the corresponding cellularity voxel
        dose_nda[i]=d

    dose_img = sitk.GetImageFromArray(dose_nda) #Convert the numpy array into an image 
    dose_img.CopyInformation(cell_map) #Copy the headers of the cellularity image onto the headers of the optimised dose image
    
    sitk.WriteImage(dose_img, out_dir+'/Dose_optimised.nii') #Save the prescribed optimised dose image as a nifti file
