# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 10:11:57 2023

@author: cbri3325
"""


#%% Import functions 

import SimpleITK as sitk
import os
import numpy as np
from ImageAnalysisFunctions import *

#%% Set Working directory

# #SBC        
# data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory
# subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
# subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#SC
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/nifti/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('P29')
subjs_name.remove('P35') #no ptv
subjs_name.remove('P39')
subjs_name.remove('P41')
subjs_name.remove('P45')
subjs_name.remove('P50')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    print('Checking maximum dose prescriptions for '+current)
    subj_dir = data_supradir+current
    subj_name = current
    
    #Read directory's path
    prsc_dir = data_supradir+current+'/RTPRESC'
    
    #Read dose prescription image
    DPresc = sitk.ReadImage(prsc_dir +'/DP_pers_CfBoost.nii')
    
    maxVal, seed = getMaxVox(DPresc)
    print('Max dose for '+current+' is '+str(maxVal))
    
    DPresc_arr = sitk.GetArrayFromImage(DPresc)
    n_vox_max = np.size(DPresc_arr[DPresc_arr>90])
    tot_vox = np.size(DPresc_arr)
    perc = (n_vox_max/tot_vox)*100
    print(str(n_vox_max)+' voxels out of '+str(tot_vox)+' ('+str(perc)+'%) have a dose prescription > 90 Gy')