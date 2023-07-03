# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 10:41:11 2023


@author: Caterina Brighi

This script allows to calculate stats on dose prescriptions/plans obtained in CTV for Chordoma patients
"""

import SimpleITK as sitk
import os
import pandas as pd
import glob
from ImageAnalysisFunctions import generate_mask

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

def Resample_image(input_image, reference_image, interpolator):
    
    '''Returns the input image resampled to the reference image space.
        interpolator can be selected amongst: sitk.sitkLinear, sitk.sitkNearestNeighbor or sitk.sitkBSpline
       Remember to write the output image into an image file after applying this function'''
       
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(reference_image)   
    resample.SetInterpolator(interpolator)
    output_image = resample.Execute(input_image)
    return output_image

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
# subjs_name=['AIRC24946_R012']

#Create an excel spreadsheet with the stats results
Results = pd.ExcelWriter(data_supradir +'Dose_stats_results.xlsx')
RBE_stats_df = pd.DataFrame(columns=['Patient_ID', 'Volume [mm3]', 'Mean dose', 'Std mean dose', 'Median dose', 'Max dose', 'Min dose'])
PHYS_stats_df = pd.DataFrame(columns=['Patient_ID', 'Volume [mm3]', 'Mean dose', 'Std mean dose', 'Median dose', 'Max dose', 'Min dose'])
Dp_stats_df = pd.DataFrame(columns=['Patient_ID', 'Volume [mm3]', 'Mean dose', 'Std mean dose', 'Median dose', 'Max dose', 'Min dose'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/' #Set path to subject directory
    
    #Read CTV, RBE dose and Dose prescription 
    
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'PTV*'):
        if (('PTV_5' in filename) or ('PTV5' in filename) or ('PTV_3' in filename) or ('PTV_low' in filename) or ('PTV_LD' in filename)):
            PTV_path = filename

    PTV = sitk.ReadImage(PTV_path) #read PTV LD
    
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'CTV*'):
        if (('CTV_5' in filename) or ('CTV5' in filename) or ('CTV_3' in filename) or ('CTV_low' in filename) or ('CTV_LD' in filename)):
            CTV_path = filename

    CTV = sitk.ReadImage(CTV_path) #read PTV LD
    CTV = generate_mask(CTV, PTV)
    # CTV = sitk.ReadImage(subj_dir+'MRI/baseline/orig/CTV_inDWI_ITK.nii')
    RBE = sitk.ReadImage(subj_dir+'RTDOSE/RBEdoseinCT.nii')
    CT = sitk.ReadImage(subj_dir+'ct.nii')
    PHYS = sitk.ReadImage(subj_dir+'RTDOSE/PHYS_TOT.nii')
    PHYS = Resample_image(PHYS, CT, sitk.sitkLinear)
    sitk.WriteImage(PHYS, subj_dir+'RTDOSE/PHYS_TOTinCT.nii')
    Dp = sitk.ReadImage(subj_dir+'Dose_painted/DP_dose_reorient.nii')
    
    #Extract stats from Phys and Dpresc in CTV
    RBE_stats = getStatsRoi(CTV,RBE)
    PHYS_stats = getStatsRoi(CTV,PHYS)
    Dp_stats = getStatsRoi(CTV,Dp)
    
    #Save stats in dataframe
    RBEdf_stats={'Patient_ID':current, 'Volume [mm3]':RBE_stats.get('Volume [mm3]'), 'Mean dose':RBE_stats.get('Mean intensity'), 'Std mean dose':RBE_stats.get('Std of Mean'), 'Median dose':RBE_stats.get('Median intensity'), 'Max dose':RBE_stats.get('Max intensity'), 'Min dose':RBE_stats.get('Min intensity')}
    RBE_stats_df = RBE_stats_df.append(RBEdf_stats, ignore_index=True)
    
    PHYSdf_stats={'Patient_ID':current, 'Volume [mm3]':PHYS_stats.get('Volume [mm3]'), 'Mean dose':PHYS_stats.get('Mean intensity'), 'Std mean dose':PHYS_stats.get('Std of Mean'), 'Median dose':PHYS_stats.get('Median intensity'), 'Max dose':PHYS_stats.get('Max intensity'), 'Min dose':PHYS_stats.get('Min intensity')}
    PHYS_stats_df = PHYS_stats_df.append(PHYSdf_stats, ignore_index=True)
    
    Dpdf_stats={'Patient_ID':current, 'Volume [mm3]':Dp_stats.get('Volume [mm3]'), 'Mean dose':Dp_stats.get('Mean intensity'), 'Std mean dose':Dp_stats.get('Std of Mean'), 'Median dose':Dp_stats.get('Median intensity'), 'Max dose':Dp_stats.get('Max intensity'), 'Min dose':Dp_stats.get('Min intensity')}
    Dp_stats_df = Dp_stats_df.append(Dpdf_stats, ignore_index=True)
    
RBE_stats_df.to_excel(Results, sheet_name='RBE stats', index=False)
PHYS_stats_df.to_excel(Results, sheet_name='PHYS DOSE stats', index=False)
Dp_stats_df.to_excel(Results, sheet_name='Dose painting plan stats', index=False)
Results.save()