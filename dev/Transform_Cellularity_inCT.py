# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 14:22:51 2023

@author: cbri3325
"""


#%% Import functions 

import SimpleITK as sitk
import os
import glob
from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import allVoxInt

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/nifti/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('P29')
subjs_name.remove('P39')
subjs_name.remove('P41')
subjs_name.remove('P45')
subjs_name.remove('P50')
#%%Create a for loop to transform cellularity from DWI to CT space for each subject sequentially

for current in subjs_name:
    
    print('Converting cellularity map into CT space for '+current)
    subj_name = current
    subj_dir = data_supradir+current
    MRI_dir = subj_dir+'/MRI/'

    studies = [ s.name for s in os.scandir(MRI_dir) if s.is_dir() ]
    studies = [x.replace('study_', '') for x in studies] #Remove word "study_" from list of studies
    studies2 = [int(x) for x in studies] #Convert list of strings in list of integers
    small_study = min(studies2) #Find the smallest study date
    study = 'study_'+str(small_study) #Define which one is the baseline (earliest) study

    if not os.path.exists(MRI_dir+study+'/micro_denoise/'):
        print('No images at baseline for '+current)
    else:

        
        if not os.path.exists(MRI_dir+study+'/micro_denoise/cellApp_CORRECTED.nii'):
            print('No cellApp image at baseline for '+current)
                   
        #Read directories' path
        micro_dir = MRI_dir+study+'/micro_denoise'
        reg_dir = MRI_dir+study+'/regFiles'
        str_dir = subj_dir+'/RTSTRUCT'
    
        #Read cellularity map
        cell_map = sitk.ReadImage(micro_dir+'/cellApp_CORRECTED.nii')
        
        #Read CTV and GTV voxels ok
        CTV_ok = sitk.ReadImage(micro_dir+'/ctv_voxel_ok_3D.nii.gz')
        GTV_ok = sitk.ReadImage(micro_dir+'/gtv_voxel_ok_3D.nii.gz')
        
        #Read CT image
        CT = sitk.ReadImage(subj_dir+'/ct.nii')
        
        #%%Read transformation files
        transform1 = sitk.ReadTransform(reg_dir+'/xf_CTonT2_1.txt')  
        transform2 = sitk.ReadTransform(reg_dir+'/xf_CTonT2_2.txt')
    
        composite_transform = sitk.CompositeTransform([transform2, transform1])
        inverse_transform = composite_transform.GetInverse()
        
        #%%Apply inverse transform to cell maps to bring them in CT space from DWI space
        cell_inCT = sitk.Resample(cell_map, CT, inverse_transform, sitk.sitkNearestNeighbor)
        sitk.WriteImage(cell_inCT, subj_dir+'/cell_density_inCT.nii')
        
        #%%Apply inverse transform to CTV ok to bring them in CT space from DWI space
        CTV_inCT = sitk.Resample(CTV_ok, CT, inverse_transform, sitk.sitkNearestNeighbor)
        sitk.WriteImage(CTV_inCT, subj_dir+'/CTV_ok_inCT.nii')
        
        #%%Apply inverse transform to CTV ok to bring them in CT space from DWI space
        GTV_inCT = sitk.Resample(GTV_ok, CT, inverse_transform, sitk.sitkNearestNeighbor)
        sitk.WriteImage(GTV_inCT, subj_dir+'/GTV_ok_inCT.nii')
        