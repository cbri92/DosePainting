# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 11:27:31 2023

@author: cbri3325
"""



#%% Import functions 

import SimpleITK as sitk
import os
import glob
from ImageAnalysisFunctions import *
from ConvertNii_ToDoseFiles import * #Use this functions only when using the dicom CT series file as reference dicom

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    print('Creating Boost plan dose prescriptions for '+current)
    subj_dir = data_supradir+current
    subj_name = current
    
    #Read directories' path
    prsc_dir = data_supradir+current+'/RTPRESC'
    rtstruct_dir = data_supradir+current+'/RTSTRUCT'
    
    #Read DP presc image
    DP_prsc = sitk.ReadImage(prsc_dir+'/DP_noBone_CfBoost.nii')
    
    #Read PTV LD
    for filename in glob.glob(rtstruct_dir+'/PTV*'):
        if (('PTV_5' in filename) or ('PTV5' in filename) or ('PTV_3' in filename) or ('PTV_low' in filename) or ('PTV_LD' in filename)):
            PTV_LD_path = filename
                        
    PTV_LD = sitk.ReadImage(PTV_LD_path) #read PTV LD
    
    #Read PTV HD
    for filename in glob.glob(rtstruct_dir+'/PTV*'):
        if (('PTV_7' in filename) or ('PTV7' in filename) or ('PTV_high' in filename) or ('PTV_HD' in filename)):
            PTV_HD_path = filename
                        
    PTV_HD = sitk.ReadImage(PTV_HD_path) #read PTV HD
    
    #Generate margins
    margin = (PTV_LD-PTV_HD)>0
    sitk.WriteImage(margin, prsc_dir+'/margins.nii')
    
    #%%Assign 74 Gy to all voxels in HD_margins, and 54 Gy to all voxels in LD_margins
    Boost_prsc_inPTV = set_mask_value(DP_prsc, PTV_HD, 74)  
    Boost_prsc_inPTV = set_mask_value(Boost_prsc_inPTV, margin, 54)
    
    sitk.WriteImage(Boost_prsc_inPTV, prsc_dir+'/Boost_presc_inPTV.nii')