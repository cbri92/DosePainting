# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 11:27:31 2023

@author: Caterina Brighi

This script generates dose prescription maps for the boost baseline plans in the CTV, receiving a uniform dose of radiation of 74 Gy.
"""

#%% Import functions 

import SimpleITK as sitk
import os
import glob
from ImageAnalysisFunctions import set_mask_value
from ConvertNii_ToDoseFiles import * #Use this functions only when using the dicom CT series file as reference dicom

#%% Set Working directory
        
data_supradir = 'path/to/pathients/data/supra/directory/' #Set working directory

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
    DP_prsc = sitk.ReadImage(prsc_dir+'/Dpainted.nii')

    #Read CTV HD
    for filename in glob.glob(rtstruct_dir+'/CTV*'):
        if (('CTV_7' in filename) or ('CTV7' in filename) or ('CTV_high' in filename) or ('CTV_HD' in filename)):
            CTV_HD_path = filename
                        
    CTV_HD = sitk.ReadImage(CTV_HD_path) #read PTV HD

    #%%Assign 74 Gy to all voxels in HD_margins, and 54 Gy to all voxels in LD_margins
    Boost_prsc_inCTV = set_mask_value(DP_prsc, CTV_HD, 74)  
    
    sitk.WriteImage(Boost_prsc_inCTV, prsc_dir+'/Boost_presc_inCTV_HD.nii')
