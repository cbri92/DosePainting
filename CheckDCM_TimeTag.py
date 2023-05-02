# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:16:10 2023

@author: cbri3325

This script checks that the scan time of the CT and RTSTRUCT dcm files match
"""

import os
import glob
import pydicom

RT_dir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/dicom/RTDATA/'

subjs_path = [ f.path for f in os.scandir(RT_dir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(RT_dir) if f.is_dir() ] #Create a list of subjects names

for current in subjs_name:
    
    print('for patient '+current)
    study_path = RT_dir+current
   # study_path = [ f.path for f in os.scandir(RT_dir+current) if f.is_dir() ][0]

    CT_dir = study_path+'/CT'
    RTSTRUCT_dir = study_path+'/RTSTRUCT'
        
    CT_dcm_file = glob.glob(CT_dir+"/*.dcm")[0] # This will get each dicom file for this patient and sequence
    CT_ds = pydicom.read_file(CT_dcm_file) # pydicom should be able to now read these Path objects, if not try updating your pydicom version :)
    print(CT_ds.StudyDate)
    
    RTSTRUCT_dcm_file = glob.glob(RTSTRUCT_dir+"/*.dcm")[0] # This will get each dicom file for this patient and sequence
    RTSTRUCT_ds = pydicom.read_file(RTSTRUCT_dcm_file) # pydicom should be able to now read these Path objects, if not try updating your pydicom version :)
    print(RTSTRUCT_ds.StudyDate)
    
    if not (CT_ds.StudyDate == RTSTRUCT_ds.StudyDate):
        print('RTSTRUCT from different CT scan for patient '+current)
    else:
        print('ok')
