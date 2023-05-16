# -*- coding: utf-8 -*-
"""
Created on Tue May  2 14:10:29 2023

@author: Caterina Brighi
"""

import numpy as np
import SimpleITK as sitk
import os
from ImageAnalysisFunctions import Resample_image
from ImageStatisticsFunctions import allVoxInt
import pandas as pd
import math

#%%Define constants

a_b = 2#Gy-1
r0 = 10**7#cm-3

#%%Functions

def find_Dlq(D_98, a_b, n):
    
    ''''Calculates Dlq from D_98, a/b and n.
    
    Inputs:   
        D_98 = dose delivered to 98% of GTV volume (float)
        a_b = alpha on beta ratio (float)
        n = number of fractions (float)
    
    Output:
        Dlq (float)'''
    
    Dlq = D_98 + D_98*(D_98/n)*(1/a_b)
    return Dlq

def find_N0(r0, roi):
    
    ''''This fucntion calculates the number of cells in a specific roi, given the cell density r0 and the volume as a sitk Image.
     
    Inputs:   
        r0 = cellular density in cm-3
        roi = volume as a sitk binary image
    
    Output:
        N0 = total number of cells in roi (float)'''
    
    voxel_volume = (roi.GetSpacing()[0]*roi.GetSpacing()[1]*roi.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
    nda = sitk.GetArrayFromImage(roi)
    volume = np.count_nonzero(nda)*voxel_volume # volume of ROI in cm-3
    N0 = r0*volume #N of cells in the volume
    
    return N0

def find_D_98(RBE, GTV):
    
    '''This function finds D_98 given a RBE dose file and a GTV volume.
    
    Inputs:        
        RBE = Relative biological dose image. sitk Image
        GTV = gross tumour volume. sitk Image
    
    Output:
        D_98 = dose delivered to 98% of GTV volume (float)'''
    
    x = allVoxInt(RBE, GTV) #This function calculates a 2D flat array of the dose for each voxel within the 3D structure
    voxInBin = math.floor(0.02*len(x)) #number of voxels in each bin as 2% of the GTV volume
    n_bins = round(len(x)/voxInBin)
    counts, bin_edges = np.histogram(x, bins=n_bins, range=(0, x.max()), normed=None, weights=None, density=False)
    histcum = 100*(1 - np.cumsum(counts)/len(x)) #cumulative histogram values: y axis
    dvh_dict = {'Dose [Gy]': bin_edges[:-1], 'Relative Volume [%]': histcum}
    DVH_df = pd.DataFrame(data=dvh_dict)    
    df_closest = DVH_df.iloc[(DVH_df['Relative Volume [%]']-98).abs().argsort()[:1]]
    Vol_98 = df_closest['Relative Volume [%]'].tolist()[0]
    D_98 = df_closest['Dose [Gy]'].tolist()[0]
    print('Relative Volume [%]: '+str(Vol_98)+ ' and Dose [Gy]: '+str(D_98))
    
    return D_98

data_dir='C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/nifti/'

Results = pd.ExcelWriter('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/DVH_info.xlsx')
dvh_df = pd.DataFrame(columns=['Patient_ID', 'N0', 'D_98% [Gy]', 'Dlq [Gy]'])

subjs_path = [f.path for f in os.scandir(data_dir) if f.is_dir()] #Create a list of the paths to the subjects directories
subjs_name = [f.name for f in os.scandir(data_dir) if f.is_dir()] #Create a list of subjects names

n_subj = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name: 
    
    subj_dir = data_dir+current    
    subj_name = current
     
    print(current) 

    # Read rbe images and generate rbe tot image
    rbe_0 = sitk.ReadImage(subj_dir + '/RTDOSE/RBE_v0.nii')
    rbe_1 = sitk.ReadImage(subj_dir + '/RTDOSE/RBE_v1.nii')
    RBE_tot = rbe_0+rbe_1
    sitk.WriteImage(RBE_tot, subj_dir + '/RTDOSE/RBE_tot.nii')
    
    # Resample RBE_tot image into ct resolution
    ct = sitk.ReadImage(subj_dir + '/ct.nii')
    RBE = Resample_image(RBE_tot, ct, sitk.sitkLinear)
    sitk.WriteImage(RBE, subj_dir + '/RTDOSE/RBE_inCT.nii')
    
    # Read GTV image  
    GTV = sitk.ReadImage(subj_dir + '/RTSTRUCT/GTV.nii.gz')
    GTV.SetOrigin(ct.GetOrigin())
    GTV.SetDirection(ct.GetDirection())
    
    # Calculate D_98
    D_98 = find_D_98(RBE,GTV)
    
    # Calculate Dlq
    Dlq = find_Dlq(D_98, a_b, 16)
    
    # Calculate N0
    N0 = find_N0(r0, GTV)
    
    #Save stats in dataframe
    dvh_dict={'Patient_ID':current, 'N0':N0, 'D_98% [Gy]':D_98, 'Dlq [Gy]':Dlq}
    dvh_df = dvh_df.append(dvh_dict, ignore_index=True)
    
#%%On local control patients

data_dir='C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/nifti_localcontrol/'

subjs_path = [f.path for f in os.scandir(data_dir) if f.is_dir()] #Create a list of the paths to the subjects directories
subjs_name = [f.name for f in os.scandir(data_dir) if f.is_dir()] #Create a list of subjects names

n_subj = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name: 
    
    subj_dir = data_dir+current    
    subj_name = current
     
    print(current) 

    # Read rbe images and generate rbe tot image
    rbe_0 = sitk.ReadImage(subj_dir + '/RTDOSE/RBE_v0.nii')
    rbe_1 = sitk.ReadImage(subj_dir + '/RTDOSE/RBE_v1.nii')
    
    if rbe_0.GetSpacing()[0] != rbe_1.GetSpacing()[0]:
        print('Different spacing')
        if rbe_0.GetSpacing()[0] > rbe_1.GetSpacing()[0]:
            print('rbe_0 is greater than rbe_1')
            rbe_1 = Resample_image(rbe_1, rbe_0, sitk.sitkLinear)
        elif rbe_0.GetSpacing()[0] < rbe_1.GetSpacing()[0]:
            print('rbe_1 is greater than rbe_0')
            rbe_0 = Resample_image(rbe_0, rbe_1, sitk.sitkLinear)
    
    RBE_tot = rbe_0+rbe_1
    sitk.WriteImage(RBE_tot, subj_dir + '/RTDOSE/RBE_tot.nii')
    
    # Resample RBE_tot image into ct resolution
    ct = sitk.ReadImage(subj_dir + '/ct.nii')
    RBE = Resample_image(RBE_tot, ct, sitk.sitkLinear)
    sitk.WriteImage(RBE, subj_dir + '/RTDOSE/RBE_inCT.nii')
    
    # Read GTV image  
    GTV = sitk.ReadImage(subj_dir + '/RTSTRUCT/GTV.nii.gz')
    GTV.SetOrigin(ct.GetOrigin())
    GTV.SetDirection(ct.GetDirection())
    
    # Calculate D_98
    D_98 = find_D_98(RBE,GTV)
    
    # Calculate Dlq
    Dlq = find_Dlq(D_98, a_b, 16)
    
    # Calculate N0
    N0 = find_N0(r0, GTV)
    
    #Save stats in dataframe
    dvh_dict={'Patient_ID':current, 'N0':N0, 'D_98% [Gy]':D_98, 'Dlq [Gy]':Dlq}
    dvh_df = dvh_df.append(dvh_dict, ignore_index=True)
    
    
    
dvh_df.to_excel(Results, sheet_name='DVH stats', index=False)
Results.save()
