# -*- coding: utf-8 -*-
"""
Created on Fri Dec 2 11:01:47 2022

@author: Caterina Brighi

This script allows to derive maps of alpha and beta from maps of LET for skull based chordoma patients treated with proton therapy.
The implementation leverages the modelling published in:
    
McNamara AL, Schuemann J, Paganetti H. A phenomenological relative biological effectiveness (RBE) model for proton therapy 
based on all published in vitro cell survival data. Phys Med Biol [Internet]. 2015 Nov 7;60(21):8399–416. 
Available from: https://iopscience.iop.org/article/10.1088/0031-9155/60/21/8399
"""

import numpy as np
import math
import SimpleITK as sitk
import os

def generate_mask(image, roi):
    
    '''Returns the masked image of the roi applied to the image.
    Remember to save the masked image file after applying this function.'''
    
    masking = sitk.MaskImageFilter() 
    mask = masking.Execute(image, roi)
    return mask

def Resample_image(input_image, reference_image):
    
    '''Returns the input image resampled to the reference image space.
       Remember to write the output image into an image file after applying this function'''
       
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(reference_image)   
    output_image = resample.Execute(input_image)
    return output_image

#Set modelling coefficients
p0 = 0.99064
p1 = 0.35605
p2 = 1.1012
p3 = -0.0038703

#Set alpha/beta photons
alpha_beta_x = 2 #Gy for chordomas
alpha_x = 0.1 #Gy-1
beta_x = alpha_x/alpha_beta_x #Gy-2

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/' #Set path to subject directory
    
    #Read CTV on CT
    CTV = sitk.ReadImage(subj_dir+'RTSTRUCT/CTV_AIRC24946.nii')
    
    #%%Read LET and PHYS dose maps 
    LET1 = sitk.ReadImage(subj_dir+'RTDOSE/LET_v0.nii')
    LET2 = sitk.ReadImage(subj_dir+'RTDOSE/LET_v1.nii')
    PHYS1 = sitk.ReadImage(subj_dir+'RTDOSE/PHYS_v0.nii')
    PHYS2 = sitk.ReadImage(subj_dir+'RTDOSE/PHYS_v1.nii')
    
    #Read CT
    CT = sitk.ReadImage(subj_dir+'ct.nii')
    
    # #Mask the LET and PHYS maps to CTV
    # LET1 = generate_mask(LET1, CTV)
    # LET2 = generate_mask(LET2, CTV)
    # PHYS1 = generate_mask(PHYS1, CTV)
    # PHYS2 = generate_mask(PHYS2, CTV)
    
    #Calculate the dose-weighted LET map
    LET = ((PHYS1*LET1)+(PHYS2*LET2))/(PHYS1+PHYS2) #From Matsumoto S, Lee SH, Imai R, Inaniwa T, Matsufuji N, Fukahori M, et al. Unresectable chondrosarcomas treated with carbon ion radiotherapy: Relationship between dose-averaged linear energy transfer and local recurrence. Anticancer Res. 2020;40(11):6429–35.
    
    #Resize dose-weighted LET to CT resolution
    LET = Resample_image(LET, CT)   
    sitk.WriteImage(LET, subj_dir+'RTDOSE/LETdoseinCT.nii')

    #%%Derive RBEmax and RBEmin using modelling from McNamara et al. 2015
    RBEmax = p0 + (p1/alpha_beta_x)*LET
    RBEmin = p2 + p3*math.sqrt(alpha_beta_x)*LET
    
    # RBEmax = generate_mask(RBEmax, CTV)
    # RBEmin = generate_mask(RBEmin, CTV)
    
    #%%Derive alpha and beta proton using estimates of RBEmax and RBEmin
    alpha = (RBEmax*alpha_x)/LET
    beta = (beta_x*(RBEmin**2))/LET
    
    #Save alpha and beta maps
    sitk.WriteImage(alpha, subj_dir+'RTDOSE/alpha.nii')
    sitk.WriteImage(alpha, subj_dir+'RTDOSE/beta.nii')
    
    #%%Apply inverse registration to map alpha and beta maps to MRI space

    transform1 = sitk.ReadTransform(subj_dir+'regFiles/xf_CTonT2.txt')
    versor1 = sitk.VersorRigid3DTransform(transform1)
    
    transform2 = sitk.ReadTransform(subj_dir+'regFiles/xf_T2onDWI.txt')
    versor2 = sitk.VersorRigid3DTransform(transform2)

    composite_transform = sitk.CompositeTransform([versor2, versor1])
    ADC = sitk.ReadImage(subj_dir+'MRI/baseline/micro/ADC_50_400_1000.nii.gz')

    alpha_reg = sitk.Resample(alpha, ADC, composite_transform, sitk.sitkNearestNeighbor)
    beta_reg = sitk.Resample(beta, ADC, composite_transform, sitk.sitkNearestNeighbor)
    CTV_reg = sitk.Resample(CTV, ADC, composite_transform, sitk.sitkNearestNeighbor)
    LET_reg = sitk.Resample(LET, ADC, composite_transform, sitk.sitkNearestNeighbor)
    CT_reg = sitk.Resample(CT, ADC, composite_transform, sitk.sitkNearestNeighbor)  
    
    sitk.WriteImage(alpha_reg, subj_dir+'RTDOSE/alpha_regDWI.nii')
    sitk.WriteImage(beta_reg, subj_dir+'RTDOSE/beta_regDWI.nii')
    sitk.WriteImage(CTV_reg, subj_dir+'MRI/baseline/micro/CTV_inDWI_ITK.nii')
    sitk.WriteImage(LET_reg, subj_dir+'RTDOSE/LETdoseinDWI.nii')
    sitk.WriteImage(CT_reg, subj_dir+'CTonDWI.nii')
    
    #%%Generate CTV in DWI excluding bones
    
    bone_mask = CT_reg>100
    bone_mask = CTV_reg*bone_mask
    CTV_noBone = CTV_reg-bone_mask
    
    sitk.WriteImage(CTV_noBone, subj_dir+'MRI/baseline/micro/CTV_inDWI_noBone.nii')
    
    #%%Mask alpha and beta maps in MRI space to CTV in DWI space

    #Read CTV on DWI
        
    alpha_DWI = generate_mask(alpha_reg, CTV_reg)
    beta_DWI = generate_mask(beta_reg, CTV_reg)
    
    alpha_DWI_noBone = generate_mask(alpha_reg, CTV_noBone)
    beta_DWI_noBone = generate_mask(beta_reg, CTV_noBone)
    
    sitk.WriteImage(alpha_DWI, subj_dir+'MRI/baseline/micro/alpha_CTVinDWI.nii')
    sitk.WriteImage(beta_DWI, subj_dir+'MRI/baseline/micro/beta_CTVinDWI.nii')
    
    sitk.WriteImage(alpha_DWI_noBone, subj_dir+'MRI/baseline/micro/alpha_CTVinDWI_noBone.nii')
    sitk.WriteImage(beta_DWI_noBone, subj_dir+'MRI/baseline/micro/beta_CTVinDWI_noBone.nii')
