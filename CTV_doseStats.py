# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 10:41:11 2023


@author: Caterina Brighi

This script allows to calculate stats on dose prescriptions/plans obtained in CTV for Chordoma patients
"""

import SimpleITK as sitk
import os
import pandas as pd

def getStatsRoi(roi, image):
    
    '''Calculates general stats from the roi applied on the image and returns them in a dictionary.
    The returned dictionary contains the follwoing parameters:
        'Volume [mm3]'
        'Mean intensity [SUV]'
        'Std of Mean [SUV]'
        'Median intensity [SUV]'
        'Max intensity [SUV]' 
        'Min intensity [SUV]'
        '''
    
    stats = sitk.LabelIntensityStatisticsImageFilter()
    stats.Execute(roi, image)
    volume = stats.GetPhysicalSize(1)
    mean = stats.GetMean(1)
    std = stats.GetStandardDeviation(1)
    median = stats.GetMedian(1)
    maxVal = stats.GetMaximum(1)
    minVal = stats.GetMinimum(1)
    return {'Volume [mm3]':volume, 'Mean intensity':mean, 'Std of Mean': std, 'Median intensity': median, 'Max intensity':maxVal, 'Min intensity':minVal}


#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#Create an excel spreadsheet with the stats results
Results = pd.ExcelWriter(data_supradir +'Dose_stats_results.xlsx')
PHYS_stats_df = pd.DataFrame(columns=['Patient_ID', 'Volume [mm3]', 'Mean dose', 'Std mean dose', 'Median dose', 'Max dose', 'Min dose'])
Dp_stats_df = pd.DataFrame(columns=['Patient_ID', 'Volume [mm3]', 'Mean dose', 'Std mean dose', 'Median dose', 'Max dose', 'Min dose'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/' #Set path to subject directory
    
    #Read CTV, Physical dose and Dose prescription on DWI
    CTV = sitk.ReadImage(subj_dir+'MRI/baseline/orig/CTV_inDWI_ITK.nii')
    PHYS = sitk.ReadImage(subj_dir+'RTDOSE/PHYS_TOTinDWI.nii')
    Dp = sitk.ReadImage(subj_dir+'MRI/baseline/Prescribed_dose/Dose_optimised_CTVnoBone_final.nii')
    
    #Extract stats from Phys and Dpresc in CTV
    PHYS_stats = getStatsRoi(CTV,PHYS)
    Dp_stats = getStatsRoi(CTV,Dp)
    
    #Save stats in dataframe
    PHYSdf_stats={'Patient_ID':current, 'Volume [mm3]':PHYS_stats.get('Volume [mm3]'), 'Mean dose':PHYS_stats.get('Mean intensity'), 'Std mean dose':PHYS_stats.get('Std of Mean'), 'Median dose':PHYS_stats.get('Median intensity'), 'Max dose':PHYS_stats.get('Max intensity'), 'Min dose':PHYS_stats.get('Min intensity')}
    PHYS_stats_df = PHYS_stats_df.append(PHYSdf_stats, ignore_index=True)
    
    Dpdf_stats={'Patient_ID':current, 'Volume [mm3]':Dp_stats.get('Volume [mm3]'), 'Mean dose':Dp_stats.get('Mean intensity'), 'Std mean dose':Dp_stats.get('Std of Mean'), 'Median dose':Dp_stats.get('Median intensity'), 'Max dose':Dp_stats.get('Max intensity'), 'Min dose':Dp_stats.get('Min intensity')}
    Dp_stats_df = Dp_stats_df.append(Dpdf_stats, ignore_index=True)
    
PHYS_stats_df.to_excel(Results, sheet_name='PHYS stats', index=False)
Dp_stats_df.to_excel(Results, sheet_name='Dose presc stats', index=False)
Results.save()