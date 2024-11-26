# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:18:20 2023

@author: cbri3325
"""


import SimpleITK as sitk
import os
import glob


def Resample_image(input_image, reference_image):
    
    '''Returns the input image resampled to the reference image space.
       Remember to write the output image into an image file after applying this function'''
       
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(reference_image)   
    output_image = resample.Execute(input_image)
    return output_image


#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SBC_Carbon/nifti_MRI/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/' #Set path to subject directory
    
    #%%Read RBE dose maps 
    
    #Read CT
    CT = sitk.ReadImage(subj_dir+'ct.nii')
    
    RBE_list = []
    for file in glob.glob(subj_dir+'RTDOSE/'+'*EFFECTIVE*'):
        RBE_list.append(file)
            
    if (not os.path.exists(subj_dir+'RTDOSE/RBEdoseinCT.nii')) and (len(RBE_list)>1):
        RBE1 = sitk.ReadImage(RBE_list[0])
        RBE1 = Resample_image(RBE1, CT) 
        RBE2 = sitk.ReadImage(RBE_list[1])
        RBE2 = Resample_image(RBE2, CT) 
        RBE = RBE1+RBE2
        sitk.WriteImage(RBE, subj_dir+'RTDOSE/RBEdoseinCT.nii')
    elif (not os.path.exists(subj_dir+'RTDOSE/RBEdoseinCT.nii')) and (len(RBE_list)==1):        
        RBE = sitk.ReadImage(RBE_list[0])
        RBE = Resample_image(RBE, CT)
        sitk.WriteImage(RBE, subj_dir+'RTDOSE/RBEdoseinCT.nii')

    
