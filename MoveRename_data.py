# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 13:20:44 2023

@author: cbri3325
"""

#%% Copy directory and its content
import shutil
import os

    
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/' #Set working directory

#%%Copy dicoms
dcm_dir = data_supradir+'dicom/'

new_dir = dcm_dir+'RTSTRUCT aggiornate/'
RT_dir = dcm_dir+'RTDATA/'

subjs_name = [ f.name for f in os.scandir(new_dir) if f.is_dir() ] #Create a list of subjects names

for subj in subjs_name:

    # #Rename CT and RTSTRUCT old folders
    # os.rename(RT_dir+subj+'/CT', RT_dir+subj+'/CT_old')
    # os.rename(RT_dir+subj+'/RTSTRUCT', RT_dir+subj+'/RTSTRUCT_old')
    
    #Copy CT and RTSTRUCT new folder in RTDATA folder
    study = [ f.name for f in os.scandir(new_dir+subj) if f.is_dir() ][0]
    shutil.copytree(new_dir+subj+'/'+study+'/CT', RT_dir+subj+'/CT')
    shutil.copytree(new_dir+subj+'/'+study+'/RTSTRUCT', RT_dir+subj+'/RTSTRUCT')


#%%Copy nifti

nii_dir = data_supradir+'nifti/'
micro_dir = data_supradir+'micro/'
microDen_dir = data_supradir+'micro_denoise/'

subjs_micro = [ f.name for f in os.scandir(micro_dir) if f.is_dir() ]

for subj in subjs_micro:
    
    studies = [ s.name for s in os.scandir(nii_dir+subj+'/MRI/') if s.is_dir() ]
    studies = [x.replace('study_', '') for x in studies] #Remove word "study_" from list of studies
    studies2 = [int(x) for x in studies] #Convert list of strings in list of integers
    small_study = min(studies2) #Find the smallest study date
    study = 'study_'+str(small_study) #Define which one is the baseline (earliest) study
    
    shutil.copytree(micro_dir+subj, nii_dir+subj+'/MRI/'+study+'/micro')
    
subjs_microDen = [ f.name for f in os.scandir(microDen_dir) if f.is_dir() ]

for subj in subjs_microDen:
    
    studies = [ s.name for s in os.scandir(nii_dir+subj+'/MRI/') if s.is_dir() ]
    studies = [x.replace('study_', '') for x in studies] #Remove word "study_" from list of studies
    studies2 = [int(x) for x in studies] #Convert list of strings in list of integers
    small_study = min(studies2) #Find the smallest study date
    study = 'study_'+str(small_study) #Define which one is the baseline (earliest) study
    
    shutil.copytree(microDen_dir+subj, nii_dir+subj+'/MRI/'+study+'/micro_denoise')