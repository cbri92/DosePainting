# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 09:50:05 2023

@author: Caterina Brighi

This script:
    1. converts boost dose prescriptions into inverse boost dose prescriptions
    2. converts inverse dose prescriptions from nifti to dicom format
"""

#%% Import functions 

import SimpleITK as sitk
import os
import glob
from ImageAnalysisFunctions import generate_mask, flip_image
from ConvertNii_ToDoseFiles import * #Use this functions only when using the dicom CT series file as reference dicom

#%% Set Working directory
        
data_supradir = 'path/to/pathients/data/supra/directory/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    print('Creating inverse dose prescriptions for '+current)
    subj_dir = data_supradir+current
    subj_name = current

    prsc_dir = data_supradir+current+'/RTPRESC'
    Boost_prsc = sitk.ReadImage(prsc_dir+'/Boost_presc_inCTV_HD.nii')
    
    #%%Generate inverse dose prescription as Dinv = Dmax - Dpresc, where Dmax = maximum dose estimated for each patient
    Inv_Boost_prsc_inCT = 81-Boost_prsc
    
    #%%Read CTV HD
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'CTV*'):
        if (('CTV_7' in filename) or ('CTV7' in filename) or ('CTV_high' in filename) or ('CTV_HD' in filename)):
            CTV_HD_path = filename
                        
    CTV_HD = sitk.ReadImage(CTV_HD_path) #read PTV HD
    
    #%%Mask inverse dose prescriptions only in voxels with valid dose prescription     
    Inv_Boost_prsc_inCTV = generate_mask(Inv_Boost_prsc_inCT, CTV_HD)
    
    #%%Write dose prescriptions and inverse dose prescriptions        
    sitk.WriteImage(Inv_Boost_prsc_inCTV, prsc_dir +'/Inv_Boost_prsc_inCTV_HD.nii')
    
    #%%Convert inverse dose prescriptions into dcm RT dose files using CT dicom series as reference dicom
    
    #First we need to flip the inverse prescriptins wrt y axis, as the dii2nii algo used originally to convert all the images flipped this wrt to original dicom images
    Inv_Boost_prsc_inCTV = flip_image(Inv_Boost_prsc_inCTV, 'y')
    ct_ref = subj_dir+'/dicom/CT/' #Path to directory of CT dicom series
    
    convert_nii_to_dicom_RTdosefile(Inv_Boost_prsc_inCTV, ct_ref, output_directory=prsc_dir, out_filename="Inv_Boost_prsc.dcm")
