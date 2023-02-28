# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:18:20 2023

@author: cbri3325
"""


import SimpleITK as sitk
import os


def Resample_image(input_image, reference_image):
    
    '''Returns the input image resampled to the reference image space.
       Remember to write the output image into an image file after applying this function'''
       
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(reference_image)   
    output_image = resample.Execute(input_image)
    return output_image


#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('AIRC24946_R052')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/' #Set path to subject directory
    
    #%%Read RBE dose maps 
    RBE = sitk.ReadImage(subj_dir+'RTDOSE/RBE_TOT.nii')
    
    #Read CT
    CT = sitk.ReadImage(subj_dir+'ct.nii')
    
    #Resize dose-weighted LET to CT resolution
    RBE = Resample_image(RBE, CT)   
    sitk.WriteImage(RBE, subj_dir+'RTDOSE/RBEdoseinCT.nii')
    
    #%%Apply inverse registration to map alpha and beta maps to MRI space

    transform1 = sitk.ReadTransform(subj_dir+'regFiles/xf_CTonT2.txt')
    versor1 = sitk.VersorRigid3DTransform(transform1)
    
    transform2 = sitk.ReadTransform(subj_dir+'regFiles/xf_T2onDWI.txt')
    versor2 = sitk.VersorRigid3DTransform(transform2)

    composite_transform = sitk.CompositeTransform([versor2, versor1])
    ADC = sitk.ReadImage(subj_dir+'MRI/baseline/micro/ADC_50_400_1000.nii.gz')

    RBE_reg = sitk.Resample(RBE, ADC, composite_transform, sitk.sitkNearestNeighbor)
    sitk.WriteImage(RBE_reg, subj_dir+'RTDOSE/RBEdoseinDWI.nii')
