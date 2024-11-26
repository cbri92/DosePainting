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
import math
import numpy as np
from ImageStatisticsFunctions import allVoxInt

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

def find_D_X(dose, percentage_volume, ROI):
    
    '''This function finds D_X given a dose file and a ROI volume.
    
    Inputs:        
        dose = dose image. sitk Image
        percentage_volume = percentage volume in which dose is calculated. int
        ROI = target volume. sitk Image
    
    Output:
        D_x = dose delivered to x% of target volume (float)'''
    
    x = allVoxInt(dose, ROI) #This function calculates a 2D flat array of the dose for each voxel within the 3D structure
    voxInBin = math.floor(0.01*len(x)) #number of voxels in each bin as 2% of the ROI volume
    n_bins = round(len(x)/voxInBin)
    counts, bin_edges = np.histogram(x, bins=n_bins, range=(0, x.max()), normed=None, weights=None, density=False)
    histcum = 100*(1 - np.cumsum(counts)/len(x)) #cumulative histogram values: y axis
    dvh_dict = {'Dose [Gy]': bin_edges[:-1], 'Relative Volume [%]': histcum}
    DVH_df = pd.DataFrame(data=dvh_dict)    
    df_closest = DVH_df.iloc[(DVH_df['Relative Volume [%]']-percentage_volume).abs().argsort()[:1]]
    Vol_X = df_closest['Relative Volume [%]'].tolist()[0]
    D_X = df_closest['Dose [Gy]'].tolist()[0]
    print('Relative Volume [%]: '+str(Vol_X)+ ' and Dose [Gy]: '+str(D_X))   
    return D_X


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
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SBC_Tutti/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#Create an excel spreadsheet with the stats results
Results = pd.ExcelWriter(data_supradir +'CTV_GTV_Dose_stats_results.xlsx')
BL_stats_df = pd.DataFrame(columns=['Patient_ID', 'GTV Volume [mm3]', 'GTV Mean dose', 'GTV Std mean dose', 'GTV Median dose', 'GTV Max dose', 'GTV Min dose', 'GTV D95%','GTV D1%','CTV Volume [mm3]', 'CTV Mean dose', 'CTV Std mean dose', 'CTV Median dose', 'CTV Max dose', 'CTV Min dose','CTV D95%','CTV D1%'])
Dp_stats_df = pd.DataFrame(columns=['Patient_ID', 'GTV Volume [mm3]', 'GTV Mean dose', 'GTV Std mean dose', 'GTV Median dose', 'GTV Max dose', 'GTV Min dose', 'GTV D95%','GTV D1%','CTV Volume [mm3]', 'CTV Mean dose', 'CTV Std mean dose', 'CTV Median dose', 'CTV Max dose', 'CTV Min dose','CTV D95%','CTV D1%'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current #Set path to subject directory
    
    if os.path.exists(subj_dir+'/PLANS/'):
    
        rtstruct_dir = subj_dir+'/RTSTRUCT'
        prsc_dir = subj_dir+'/RTPRESC'
        
        #Read CTV HD
        for filename in glob.glob(rtstruct_dir+'/CTV*'):
            if (('CTV_7' in filename) or ('CTV_high' in filename) or ('CTV_HD' in filename)):
                CTV_HD_path = filename
                            
        CTV = sitk.ReadImage(CTV_HD_path) #read PTV LD
           
        #Read GTV
        GTV = sitk.ReadImage(rtstruct_dir+'/GTV.nii')
      
        print('Reading DOSE images for '+current) 
    
        DP_dose_path = subj_dir+'/PLANS/Dose_painting_4B/DP_4B_reorient.nii'   
        DP_dose = sitk.ReadImage(DP_dose_path) #Read DP dose
        
        BL_dose_path = subj_dir+'/PLANS/Baseline_4B/Baseline_4B.nii'
        BL_dose = sitk.ReadImage(BL_dose_path) #Read BL dose
    
        #Extract stats from BL and DP plan in GTV
        BL_stats_GTV = getStatsRoi(GTV,BL_dose)
        Dp_stats_GTV = getStatsRoi(GTV,DP_dose)
        
        #Extract stats from BL and DP plan in CTV
        BL_stats_CTV = getStatsRoi(CTV,BL_dose)
        Dp_stats_CTV = getStatsRoi(CTV,DP_dose)
        
        #Extract D95% from BL and DP plan in GTV
        BL_D95_GTV = find_D_X(BL_dose, 95, GTV)
        Dp_D95_GTV = find_D_X(DP_dose, 95, GTV)
        
        #Extract D95% from BL and DP plan in CTV
        BL_D95_CTV = find_D_X(BL_dose, 95, CTV)
        Dp_D95_CTV = find_D_X(DP_dose, 95, CTV)
        
        #Extract D1% from BL and DP plan in GTV
        BL_D1_GTV = find_D_X(BL_dose, 1, GTV)
        Dp_D1_GTV = find_D_X(DP_dose, 1, GTV)
        
        #Extract D1% from BL and DP plan in CTV
        BL_D1_CTV = find_D_X(BL_dose, 1, CTV)
        Dp_D1_CTV = find_D_X(DP_dose, 1, CTV)
        
        #Save stats in dataframe
        BLdf_stats={'Patient_ID':current, 'GTV Volume [mm3]':BL_stats_GTV.get('Volume [mm3]'), 'GTV Mean dose':BL_stats_GTV.get('Mean intensity'), 'GTV Std mean dose':BL_stats_GTV.get('Std of Mean'), 'GTV Median dose':BL_stats_GTV.get('Median intensity'), 'GTV Max dose':BL_stats_GTV.get('Max intensity'), 'GTV Min dose':BL_stats_GTV.get('Min intensity'),'GTV D95%':BL_D95_GTV, 'GTV D1%':BL_D1_GTV, 'CTV Volume [mm3]':BL_stats_CTV.get('Volume [mm3]'), 'CTV Mean dose':BL_stats_CTV.get('Mean intensity'), 'CTV Std mean dose':BL_stats_CTV.get('Std of Mean'), 'CTV Median dose':BL_stats_CTV.get('Median intensity'), 'CTV Max dose':BL_stats_CTV.get('Max intensity'), 'CTV Min dose':BL_stats_CTV.get('Min intensity'),'CTV D95%':BL_D95_CTV, 'CTV D1%':BL_D1_CTV}
        BL_stats_df = BL_stats_df.append(BLdf_stats, ignore_index=True)
          
        Dpdf_stats={'Patient_ID':current, 'GTV Volume [mm3]':Dp_stats_GTV.get('Volume [mm3]'), 'GTV Mean dose':Dp_stats_GTV.get('Mean intensity'), 'GTV Std mean dose':Dp_stats_GTV.get('Std of Mean'), 'GTV Median dose':Dp_stats_GTV.get('Median intensity'), 'GTV Max dose':Dp_stats_GTV.get('Max intensity'), 'GTV Min dose':Dp_stats_GTV.get('Min intensity'),'GTV D95%':Dp_D95_GTV, 'GTV D1%':Dp_D1_GTV, 'CTV Volume [mm3]':Dp_stats_CTV.get('Volume [mm3]'), 'CTV Mean dose':Dp_stats_CTV.get('Mean intensity'), 'CTV Std mean dose':Dp_stats_CTV.get('Std of Mean'), 'CTV Median dose':Dp_stats_CTV.get('Median intensity'), 'CTV Max dose':Dp_stats_CTV.get('Max intensity'), 'CTV Min dose':Dp_stats_CTV.get('Min intensity'),'CTV D95%':Dp_D95_CTV, 'CTV D1%':Dp_D1_CTV}
        Dp_stats_df = Dp_stats_df.append(Dpdf_stats, ignore_index=True)
        
BL_stats_df.to_excel(Results, sheet_name='Baseline 4B stats', index=False)
Dp_stats_df.to_excel(Results, sheet_name='Dose painting 4B stats', index=False)
Results.save()