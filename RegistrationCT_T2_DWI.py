# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 10:03:05 2023

@author: Caterina Brighi

This script register CT data to DWI data passing through T2 data in chordoma patients.
"""

#%% Import functions 

import SimpleITK as sitk
import glob
import os
import pandas as pd
from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import *


#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current
    if not os.path.exists(data_supradir+current+'/MRI/baseline/Prescribed_dose'):#if it does not already exist, create a directory where the optimized dose will be saved
        os.mkdir(data_supradir+current+'/MRI/baseline/Prescribed_dose')
    out_dir = data_supradir+current+'/MRI/baseline/Prescribed_dose'
    