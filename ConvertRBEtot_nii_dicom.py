# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:24:10 2023

@author: cbri3325
"""
import SimpleITK as sitk
from ConvertNii_ToDoseFiles import * #Use this functions only when using the dicom CT series file as reference dicom
from ImageAnalysisFunctions import Resample_image, flip_image


#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    print('Converting RBEtot for '+current)
    subj_dir = data_supradir+current
    subj_name = current
    
    ct_img = sitk.ReadImage(subj_dir+'/ct.nii')
    RBE_img = sitk.ReadImage(subj_dir+'/RTDOSE/RBE_TOT.nii')
    RBE_img = Resample_image(RBE_img, ct_img, sitk.sitkLinear)
    #First we need to flip the inverse prescriptins wrt y axis, as the dii2nii algo used originally to convert all the images flipped this wrt to original
    RBE_img = flip_image(RBE_img, 'y')
    
    ct_ref = subj_dir+'/CT/' #Path to directory of CT dicom series
    convert_nii_to_dicom_RTdosefile(RBE_img, ct_ref, output_directory=subj_dir+'/RTDOSE/', out_filename="RBE_TOT.dcm")
