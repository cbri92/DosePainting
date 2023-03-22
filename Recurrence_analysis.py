# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 04:06:45 2023

@author: Caterina Brighi
"""

#%% Import functions 

import matplotlib.pyplot as plt
import matplotlib as mpl
import SimpleITK as sitk
import numpy as np
import pandas as pd
import datetime
import os
import glob
import gzip
import shutil
import xlsxwriter
import time
from scipy.stats.stats import pearsonr
import dicom2nifti
from multiprocessing.pool import ThreadPool
from functools import partial
from ImageAnalysisFunctions import *

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

#Generate excel file containing voxelwise results
Results = pd.ExcelWriter(data_supradir +'Recurrence_analysis.xlsx')

#Create an empty dataframe to populate as going through the loop
df = pd.DataFrame(columns=['Subject_ID', '% Vol of recurrence in PTV', 'Mean std dose in targeted recurrence volume [Gy]', 'Max std dose in targeted recurrence volume [Gy]', 'Mean DP dose in targeted recurrence volume [Gy]', 'Max DP dose in targeted recurrence volume [Gy]', '% vol of recurrence within PTV_HD receiving 95% of prescribed HD', '% vol of recurrence within PTV_LD receiving 95% of prescribed LD'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    print('Performing recurrence dose analysis for '+current)
    subj_dir = data_supradir+current
    subj_name = current
    
    #Set path to directories
    dose_dir = subj_dir+'/RTDOSE'
    presc_dir = subj_dir+'/RTPRESC'
    struct_dir = subj_dir+'/RTSTRUCT'
    
    #Read std dose, DP presc, PTV and recurrence contours
    Std_dose = sitk.ReadImage(dose_dir+'/RBE_TOT.nii') #Read std dose image
    DPresc = sitk.ReadImage(presc_dir+'/DP_noBone_CfBoost.nii') #Read dose prescription image   
    # Rec = sitk.ReadImage(struct_dir+'/rec_per_ricerca.nii') #Read recurrence contour
    
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'PTV*'):
        if (('PTV_5' in filename) or ('PTV5' in filename) or ('PTV_3' in filename) or ('PTV_low' in filename) or ('PTV_LD' in filename)):
            PTV_LD_path = filename

    PTV_LD = sitk.ReadImage(PTV_LD_path) #read PTV LD
    
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'PTV*'):
        if (('PTV_7' in filename) or ('PTV7' in filename) or ('PTV_high' in filename) or ('PTV_HD' in filename)):
            PTV_HD_path = filename
            
    PTV_HD = sitk.ReadImage(PTV_HD_path)
    
    if current=='AIRC24946_R012' or current=='AIRC24946_R048':
        Rec = sitk.ReadImage(struct_dir+'/Recurrence_fixOrigDir.nii') #Read recurrence contour
    else:
        #Fix Origin and Direction of recurrence contour
        Rec = sitk.ReadImage(struct_dir+'/rec_per_ricerca.nii') #Read recurrence contour
        Rec = flip_image(Rec, 'y')
        Rec.SetOrigin(PTV_LD.GetOrigin())
        Rec.SetDirection(PTV_LD.GetDirection())
        sitk.WriteImage(Rec, struct_dir+'/Recurrence_fixOrigDir.nii')
    
    #Resample Std_dose image to same resolution as dose prescription
    Std_dose = Resample_image(Std_dose, DPresc, sitk.sitkLinear)
    # sitk.WriteImage(Std_dose, dose_dir+'/PHYS_TOT_inCT.nii')
    
    #%%Calculate stats for contours
    
    #Calculate % Volume of recurrence in PTV
    Tar_rec = (PTV_LD+Rec)>1 #Volume of recurrence overlapping with PTV
    sitk.WriteImage(Tar_rec, struct_dir+'/Targeted_rec.nii')
    
    Vol_targeted_rec = round(getVolumeRoi(Tar_rec, Std_dose))
    Vol_Rec = round(getVolumeRoi(Rec, Std_dose))
    Perc_targeted_rec = (Vol_targeted_rec/Vol_Rec)*100

    #Calculate mean dose in targeted recurrence volume
    Std_dose_trec_mean = getMeanRoi(Tar_rec, Std_dose)
    DPresc_trec_mean = getMeanRoi(Tar_rec, DPresc)
    
    #Calculate max dose in targeted recurrence volume
    Std_dose_trec_max = getMaxRoi(Tar_rec, Std_dose)
    DPresc_trec_max = getMaxRoi(Tar_rec, DPresc)
    
    #Calculate 95% isodose contours
    SDose_PTV_HD = generate_mask(Std_dose, PTV_HD)   
    Iso95_HD = SDose_PTV_HD>(74*0.95)
    sitk.WriteImage(Iso95_HD, struct_dir+'/Iso95_HD.nii')
    
    SDose_PTV_LD = generate_mask(Std_dose, PTV_LD)
    Iso95_LD = SDose_PTV_LD>(54*0.95)
    sitk.WriteImage(Iso95_LD, struct_dir+'/Iso95_LD.nii')
    
    #Calculate margins between Iso95 HD and LD conturs
    Iso95_HDLD_margins = (Iso95_HD+Iso95_LD)<2
    sitk.WriteImage(Iso95_HDLD_margins, struct_dir+'/Iso95_HDLD_margins.nii')
    
    #Calculate % volume of recurrence in Iso95 HD contour
    Rec_inIso95HD = (Rec+Iso95_HD)
    sitk.WriteImage(Rec_inIso95HD, struct_dir+'/Rec_inIso95HD.nii')
    Vol_Rec95HD = round(getVolumeRoi(Rec_inIso95HD, Std_dose))
    Perc_Rec95HD = (Vol_Rec95HD/Vol_Rec)*100
    
    #Calculate % volume of recurrence in Iso95 HDLD margins
    Rec_inIso95margins = (Rec+Iso95_HDLD_margins)
    sitk.WriteImage(Rec_inIso95margins, struct_dir+'/Rec_inIso95margins.nii')
    Vol_Rec95margins = round(getVolumeRoi(Rec_inIso95margins, Std_dose))
    Perc_Rec95margins = (Vol_Rec95margins/Vol_Rec)*100

    #%%Save statistics of various VOIs            
    sheet_df = {'Subject_ID': subj_name, '% Vol of recurrence in PTV':Perc_targeted_rec, 'Mean std dose in targeted recurrence volume [Gy]':Std_dose_trec_mean, 'Max std dose in targeted recurrence volume [Gy]':Std_dose_trec_max, 'Mean DP dose in targeted recurrence volume [Gy]':DPresc_trec_mean, 'Max DP dose in targeted recurrence volume [Gy]':DPresc_trec_max, '% vol of recurrence within PTV_HD receiving 95% of prescribed HD':Perc_Rec95HD, '% vol of recurrence within PTV_LD receiving 95% of prescribed LD':Perc_Rec95margins}
    df1 = pd.DataFrame(sheet_df, index=[0], columns=['Subject_ID', '% Vol of recurrence in PTV', 'Mean std dose in targeted recurrence volume [Gy]', 'Max std dose in targeted recurrence volume [Gy]', 'Mean DP dose in targeted recurrence volume [Gy]', 'Max DP dose in targeted recurrence volume [Gy]', '% vol of recurrence within PTV_HD receiving 95% of prescribed HD', '% vol of recurrence within PTV_LD receiving 95% of prescribed LD']) 
    df = df.append(df1, ignore_index=True)
    
#%%Sava data to excel spreadsheet    
df.to_excel(Results, sheet_name='Sheet1', index=False)
Results.save()