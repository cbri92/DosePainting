# -*- coding: utf-8 -*-
"""
Created on Wed May  3 17:13:55 2023

@author: cbri3325
"""
import os
import glob
import shutil
import SimpleITK as sitk
import sys
sys.path.append("C:/Users/cbri3325/Anaconda3/lib/site-packages/platipy/__init__.py") # Path containing PlatiPy library
from platipy.dicom.io.rtdose_to_nifti import convert_rtdose

data_dir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/'

pat_dir = data_dir+'AIRC24946_R012/'
DP_dir = pat_dir+'DOSEPAINTED/'
DP_nii = pat_dir+'DP_nii/'
convert_rtdose(DP_dir+'RD1.2.752.243.1.1.20230503152950410.3300.66868.dcm', False, DP_nii+'/DP.nii')
