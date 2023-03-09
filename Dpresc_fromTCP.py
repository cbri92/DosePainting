# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:31:51 2022

@author: Caterina Brighi

This script is used to take information about cellularity, alpha and beta and use them
to derive patient-specific dose prescriptions by means of a TCP model (poissonian, LQ).
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

#%%Set alpha and beta particle-cell line specific values

#Set alpha/beta photons
alpha_beta_x = 2.4 #Gy for chordomas
alpha_x = 0.1 #Gy-1
beta_x = alpha_x/alpha_beta_x #Gy-2

#Values RBEmax and RBEmin for proton on chordoma cells taken from paper Paganetti, International Journal of Radiation Oncology*Biology*Physics, 2022, 112(1), 222-236.

RBEmax = 1.59
RBEmin = 1.18

a = (RBEmax*alpha_x)
b = (beta_x*(RBEmin**2))
#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('AIRC24946_R052')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/MRI/baseline/micro/' #Set path to subject directory where cell density map image and GTV image are stored
    if not os.path.exists(data_supradir+current+'/MRI/baseline/Prescribed_dose'):#if it does not already exist, create a directory where the optimized dose will be saved
        os.mkdir(data_supradir+current+'/MRI/baseline/Prescribed_dose')
    out_dir = data_supradir+current+'/MRI/baseline/Prescribed_dose'
    
    #Unzip cellApp_CORRECTED.nii.gz file and delete original gzipped file
    for filename in glob.glob(data_supradir+current+'/MRI/baseline/orig/' +'cellApp_CORRECTED.nii.gz'):
        gunzip_shutil(filename, filename[:-3])
     
    print('Optimizing the dose prescription for '+current) 
              
    cell_map = sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/cellApp_CORRECTED.nii') #Read cellular density map. Units: um-3
    cell_map = cell_map*(10**8) #Convert cell map to Units: 10^4 cm-3
    # cell_map = cell_map*(10**12) #Convert cell map to Units: cm-3
    
    cell_dirs = [data_supradir+current+'/MRI/baseline/orig/',data_supradir+current+'/MRI/baseline/noBone/']
    cell_types = ['orig', 'noBone']
    
    for cell_dir, cell_type in zip(cell_dirs,cell_types):

        if cell_dir == data_supradir+current+'/MRI/baseline/orig/':
            # alpha_map = sitk.ReadImage(subj_dir+'alpha_CTVinDWI.nii') 
            # beta_map = sitk.ReadImage(subj_dir+'beta_CTVinDWI.nii') 
        
            CTV_full = sitk.ReadImage(cell_dir+'CTV_inDWI_ITK.nii') #Read CTV
            # CTV.SetOrigin(cell_map.GetOrigin())
            
        elif cell_dir == data_supradir+current+'/MRI/baseline/noBone/':
            # alpha_map = sitk.ReadImage(subj_dir+'alpha_CTVinDWI_noBone.nii') 
            # beta_map = sitk.ReadImage(subj_dir+'beta_CTVinDWI_noBone.nii')
        
            CTV_full = sitk.ReadImage(cell_dir+'CTV_inDWI_noBone.nii') #Read CTV
            
            cell_map=sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/cellApp_CORRECTED.nii')
            cell_map=generate_mask(cell_map, CTV_full)
            cell_map = cell_map*(10**8)
            sitk.WriteImage(cell_map, cell_dir+'cellApp_CORRECTED.nii')
            
            ctv_ok=sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/ctv_voxel_ok_3D.nii.gz', sitk.sitkUInt8)
            ctv_ok=(CTV_full+ctv_ok)>1
            sitk.WriteImage(ctv_ok, cell_dir+'ctv_voxel_ok_3D.nii.gz')
                       
        CTV = sitk.ReadImage(cell_dir+'ctv_voxel_ok_3D.nii.gz')
        CTV = sitk.Cast(CTV, sitk.sitkUInt8)
        CTV_outliers = CTV_full-CTV
        sitk.WriteImage(CTV_outliers, cell_dir+'CTV_outliers.nii')
        
        cell_map = generate_mask(cell_map, CTV)
        cell_map_nda = sitk.GetArrayFromImage(cell_map) #Convert cell map image into a numpy array
        # alpha_map_nda = sitk.GetArrayFromImage(alpha_map) #Convert alpha map image into a numpy array
        # beta_map_nda = sitk.GetArrayFromImage(beta_map) #Convert beta map image into a numpy array    
        
        #Find indexes of voxels that have values of cellularity
        result = np.where(cell_map_nda>0)
        listOfCoordinates= list(zip(result[0], result[1], result[2])) #Save indexes in a list of coordinates
        # for cord in listOfCoordinates:
        #     print(cord)    
    
        orig_shape = sitk.GetArrayFromImage(cell_map).shape #Save the shape of the original cell map image
        voxel_volume = (cell_map.GetSpacing()[0]*cell_map.GetSpacing()[1]*cell_map.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
        
        #Define input parameters for TCP modelling: number of cells per voxel (N0) and alpha and beta values per voxel
        CTV_final = CTV
        # CTV_final = (cell_map>0) & (alpha_map>0) #Create a mask where both alpha and cellularity values are >0
        
        cell_dens = allVoxInt(cell_map, CTV_final) #Define array of values with cellular density>0
        N0 = cell_dens*voxel_volume #Multiply the values of cellular density per voxel by the volume of one voxel to obtain the number of cells per voxel
        
        alpha=np.full(len(N0), a, dtype=np.float32)
        beta=np.full(len(N0), b, dtype=np.float32)
        
        # alpha = allVoxInt(alpha_map, CTV_final) #Define array of values with alpha
        # beta = allVoxInt(beta_map, CTV_final) #Define array of values with beta
        
        # alpha = np.random.normal(0.30, 0.05, size=len(N0)) #Define array of values of alpha (radiosensitivity parameter)
        # beta = np.random.normal(0.15, 0.01, size=len(N0)) #Define array of values of beta (radiosensitivity parameter)
        
        n=37 #Number of fractions
    
    #%%Optimization
        
        # for graphing below only
        x_vals = np.arange(0,100)
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
        plt.savefig(out_dir+'/TCP_dose_plot_'+cell_type+'.jpg')
        plt.show()
    
        #Plot the probability density function vs dose
        plt.plot(x_vals[:-1],np.diff(TCP))
        plt.plot([MAX_DOSE,MAX_DOSE],[0,max(np.diff(TCP))])
        plt.plot([MIN_DOSE,MIN_DOSE],[0,max(np.diff(TCP))])
        plt.legend(['PDF','Max dose value: '+str(MAX_DOSE), 'Min dose value: '+str(MIN_DOSE)], loc ="upper left")
        plt.xlabel("Dose (Gy)")
        plt.ylabel("Probability density")
        plt.savefig(out_dir+'/PDF_dose_plot_'+cell_type+'.jpg')
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
        TCP_heter =np.prod( np.exp( - N0 * np.exp(-alpha*Dopt-(beta*(Dopt**2))/n))) #TCP obtained with optimised heterogeneous dose distribution
        TCP_boost = np.prod( np.exp( - N0 * np.exp(-alpha*MAX_DOSE-(beta*(MAX_DOSE**2))/n))) # TCP with optimised homogeneous max dose boost to entire GTV
        print('TCP optimised heterogeneous dose plan: ',TCP_heter)
        print('TCP optimised homogeneous max dose boost plan: ',TCP_boost)
        
        plt.hist(Dopt, bins=20)
        plt.xlabel("Dose (Gy)")
        plt.ylabel("Counts")
        plt.savefig(out_dir+'/Dose_histo_'+cell_type+'.jpg')
        plt.close()
        
    #%% Saving the optimised dose values into a prescription dose image
        
        dose_nda = np.ndarray(shape=orig_shape) #Generate a 3D array of the same shape of the cell density map
        
        for d,i in zip(list(Dopt),listOfCoordinates): #Assign the calculated optimised dose to each voxel, based on the index of the corresponding cellularity voxel
            dose_nda[i]=d
    
        dose_img = sitk.GetImageFromArray(dose_nda) #Convert the numpy array into an image 
        dose_img.CopyInformation(cell_map) #Copy the headers of the cellularity image onto the headers of the optimised dose image
        
        sitk.WriteImage(dose_img, out_dir+'/Dose_optimised_CTV'+cell_type+'.nii') #Save the prescribed optimised dose image as a nifti file 
        
    #%%Set values within CTV with cellularity untrusted and bone pixels to receive as mean dose
        
        CTV_complete = sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/CTV_inDWI_ITK.nii')
        CTV_remaining = (CTV_complete>0) & (dose_img==0)
        mean_dose=(MAX_DOSE+MIN_DOSE)/2
        dose_img_final = set_mask_value(dose_img, CTV_remaining, mean_dose)
        sitk.WriteImage(dose_img_final, out_dir+'/Dose_optimised_CTV'+cell_type+'_final.nii')
        
        
        
