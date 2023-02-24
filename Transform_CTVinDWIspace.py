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
import glob

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('AIRC24946_R052')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/' #Set path to subject directory
    
    #Read CTV on CT
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'CTV*'):
        if (('CTV_5' in filename) or ('CTV5' in filename) or ('CTV_3' in filename) or ('CTV_low' in filename) or ('CTV_LD' in filename)):
            CTV_path = filename

    CTV = sitk.ReadImage(CTV_path) #read CTV
   
    #Read CT
    CT = sitk.ReadImage(subj_dir+'ct.nii')
    
    #%%Apply inverse registration to map alpha and beta maps to MRI space

    transform1 = sitk.ReadTransform(subj_dir+'regFiles/xf_CTonT2.txt')
    versor1 = sitk.VersorRigid3DTransform(transform1)
    
    transform2 = sitk.ReadTransform(subj_dir+'regFiles/xf_T2onDWI.txt')
    versor2 = sitk.VersorRigid3DTransform(transform2)

    composite_transform = sitk.CompositeTransform([versor2, versor1])
    ADC = sitk.ReadImage(subj_dir+'MRI/baseline/micro/ADC_50_400_1000.nii.gz')

    CTV_reg = sitk.Resample(CTV, ADC, composite_transform, sitk.sitkNearestNeighbor)
    CT_reg = sitk.Resample(CT, ADC, composite_transform, sitk.sitkNearestNeighbor)  
    
    sitk.WriteImage(CTV_reg, subj_dir+'MRI/baseline/micro/CTV_inDWI_ITK.nii')
    # sitk.WriteImage(CT_reg, subj_dir+'CTonDWI.nii')
    
    #%%Generate CTV in DWI excluding bones
    
    bone_mask = CT_reg>100
    bone_mask = CTV_reg*bone_mask
    CTV_noBone = CTV_reg-bone_mask   
    sitk.WriteImage(CTV_noBone, subj_dir+'MRI/baseline/micro/CTV_inDWI_noBone.nii')
    