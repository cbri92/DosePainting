# -*- coding: utf-8 -*-
"""
Created on Wed May  3 17:13:55 2023

@author: Caterina Brighi

This script converts dose painting and baseline dicom dose plans from dicom RT format to nifti format.
"""

import os
import pydicom
import sys
sys.path.append("C:/Users/cbri3325/Anaconda3/lib/site-packages/platipy/__init__.py") # Path containing PlatiPy library
from platipy.dicom.io.rtdose_to_nifti import convert_rtdose

#%%Set path to working directory

data_supradir = 'path/to/patients/data/supra/directory/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

n_subj = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name: 
    
    subj_dir = data_dir+current    
    subj_name = current
    
    if os.path.exists(subj_dir+'/PLANS/'):
        
        print('Converting DP 4B dose to nii for',current) 
        
        DP_dcm_dir = subj_dir+'/PLANS/DP 4B/'
        if not os.path.exists(subj_dir+'/PLANS/Dose_painting_4B/'):#if it does not already exist, create a directory where the optimized dose will be saved
            os.mkdir(subj_dir+'/PLANS/Dose_painting_4B/')
        DP_nii_dir = subj_dir+'/PLANS/Dose_painting_4B/'
        
        dcm_files = [f.path for f in os.scandir(DP_dcm_dir) if f.is_file()]
        n_files = len(dcm_files)
        for file,n in zip(dcm_files, range(n_files)):
            ds = pydicom.dcmread(file)
            if ds.Modality == 'RTDOSE':
                # print(ds.Modality)
                # print(ds.DoseType)
                # print(ds.DoseSummationType)
                # print(str(n+1))
                out_name = str(ds.SeriesDescription)+'_'+str(ds.Modality)+'_'+str(ds.DoseType)+'_'+str(ds.DoseSummationType)+'_'+str(n+1)+'.nii'
                convert_rtdose(file, False, DP_nii_dir+out_name)
    
        print('Converting Baseline 4B dose to nii for',current) 
        
        BL_dcm_dir = subj_dir+'/PLANS/Baseline/'
        if not os.path.exists(subj_dir+'/PLANS/Baseline_4B/'):#if it does not already exist, create a directory where the optimized dose will be saved
            os.mkdir(subj_dir+'/PLANS/Baseline_4B/')
        BL_nii_dir = subj_dir+'/PLANS/Baseline_4B/'
        
        dcm_files = [f.path for f in os.scandir(BL_dcm_dir) if f.is_file()]
        n_files = len(dcm_files)
        for file,n in zip(dcm_files, range(n_files)):
            ds = pydicom.dcmread(file)
            if ds.Modality == 'RTDOSE':
                # print(ds.Modality)
                # print(ds.DoseType)
                # print(ds.DoseSummationType)
                # print(str(n+1))
                out_name = str(ds.SeriesDescription)+'_'+str(ds.Modality)+'_'+str(ds.DoseType)+'_'+str(ds.DoseSummationType)+'_'+str(n+1)+'.nii'
                convert_rtdose(file, False, BL_nii_dir+out_name)
