# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 14:33:01 2021

@author: Caterina Brighi

This script derives DVH data for RT targets and organs at risk from the nifti dose and RT structures files for both the dose painting and the baseline plans.
"""

from ImageStatisticsFunctions import allVoxInt
import SimpleITK as sitk
import os.path
import sys
import numpy as np
import pandas as pd
import warnings
import shutil

#%%Set path to patients directory

data_supradir = 'path/to/patients/data/supra/directory/' #Set working directory

if not os.path.isdir(data_supradir):
    warnings.warn('invalid data_supradir supplied; quitting', stacklevel=2)
    sys.exit(1)

subjs_path = [f.path for f in os.scandir(data_supradir) if f.is_dir()] #Create a list of the paths to the subjects directories
subjs_name = [f.name for f in os.scandir(data_supradir) if f.is_dir()] #Create a list of subjects names

#Define OARs DVH files names
tronco = ["Tronco.xlsx" ,"tronco.xlsx" , "tronco_encefalico.xlsx" , "Tronco_enc.xlsx"]
chiasma = ["Chiasma.xlsx", "chiasma.xlsx"]
nervo_ottico_sx = ["Nervo_ottico_sx.xlsx", "nervo_ottico_sx.xlsx", "nervo_ott_sx.xlsx", "Nervo_ott_sx.xlsx", "n_ottico_sx.xlsx"]
nervo_ottico_dx = ["Nervo_ottico_dx.xlsx", "nervo_ottico_dx.xlsx", "nervo_ott_dx.xlsx", "Nervo_ott_dx.xlsx","n_ottico_dx.xlsx"]
coclea_sx = ["Coclea_SX.xlsx", "coclea_sx.xlsx","Coclea_sx.xlsx", "coclea_sin.xlsx"]
coclea_dx = ["Coclea_DX.xlsx", "coclea_dx.xlsx", "Coclea_dx.xlsx"]
temp_lobe_sx = ["LobeTemporal_L.xlsx", "lobo_temp_sx.xlsx"]
temp_lobe_dx = ["LobeTemporal_R.xlsx", "lobo_temp_dx.xlsx"]
carotide_sx = ["Carotide_sx.xlsx", "carotide_sx.xlsx", "carot_sx.xlsx", "carotide_sin.xlsx", "carotide_int_sx.xlsx"]
carotide_dx = ["Carotide_dx.xlsx", "carotide_dx.xlsx", "carot_dx.xlsx", "carotide_int_dx.xlsx"]
targets = ['GTV.xlsx', 'CTV_AIRC24946.xlsx']

structures_orig = tronco+chiasma+nervo_ottico_sx+nervo_ottico_dx+coclea_sx+coclea_dx+temp_lobe_sx+temp_lobe_dx+carotide_sx+carotide_dx+targets

n_subj = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    subj_dir = data_supradir+current
    
    if os.path.exists(subj_dir+'/PLANS/'):
        
        subj_name = current
        
        print('Starting DVH generation for ' + current)
        
        #%% Set path to dose and structure folders
        
        dose_dir = subj_dir + '/PLANS'
        struct_dir = subj_dir + '/RTSTRUCT'
        
        #%%Define path to dose file
        
        DP_dose = sitk.ReadImage(dose_dir+'/Dose_painting_4B/DP_4B_reorient.nii')
        BL_dose = sitk.ReadImage(dose_dir+'/Baseline_4B/Baseline_4B.nii')
        
        #%%Generate DVH folder
        if os.path.isdir(subj_dir + '/DVH_files2'):
            warnings.warn(f'emptying subj_dir directory {subj_dir}', stacklevel=2)
            shutil.rmtree(subj_dir + '/DVH_files2')
        
        os.mkdir(subj_dir + '/DVH_files')
        os.mkdir(subj_dir + '/DVH_files/Dose_painting_4B')
        os.mkdir(subj_dir + '/DVH_files/Baseline_4B')
        DP_DVH_dir = subj_dir + '/DVH_files/Dose_painting_4B'
        BL_DVH_dir = subj_dir + '/DVH_files/Baseline_4B'
        
        #%%Generate DVH excel file for each structure
        
        for f in os.scandir(struct_dir):
            file_write_name = f.name.replace(' ','_') # replcaes spaces with _
            [file_write_name, dum] = os.path.splitext(file_write_name)  # removes old extension
            file_write_name = file_write_name + '.xlsx'  # adds new extension
    
            struct = sitk.ReadImage(f.path)
            struct = struct > 0
            
            if file_write_name in structures_orig:
                
                # if file_write_name in targets:
                #     Max_dose = 74
                # elif file_write_name in tronco:
                #     Max_dose = 61
                # elif file_write_name in chiasma:
                #     Max_dose = 54
                # elif file_write_name in nervo_ottico_sx+nervo_ottico_dx:
                #     Max_dose = 50
                # elif file_write_name in coclea_sx+coclea_dx:
                #     Max_dose = 45
                # elif file_write_name in temp_lobe_sx+temp_lobe_dx:
                #     Max_dose = 71
                # elif file_write_name in carotide_sx+carotide_dx:
                #     Max_dose = 75.5

                # DP_dose = (DP_dose/Max_dose)*100
                # BL_dose = (BL_dose/Max_dose)*100
                
                #Write DVH for DP plan
                DP_x = allVoxInt(DP_dose, struct) #This function calculates a 2D flat array of the dose for each voxel within the 3D structure
                DP_counts, DP_bin_edges = np.histogram(DP_x, bins=81, range=(0, 81), weights=None, density=False)
                DP_histcum = 100*(1 - np.cumsum(DP_counts)/len(DP_x)) #cumulative histogram values: y axis
                
                #Save DP_DVH to excel file
                DP_Histo_results = pd.ExcelWriter(DP_DVH_dir +'/' + file_write_name)
                DP_df_dict = {'Dose [Gy]': DP_bin_edges[:-1], 'Relative Volume [%]': DP_histcum}
                DP_Histo_df = pd.DataFrame(data=DP_df_dict)
                DP_Histo_df.to_excel(DP_Histo_results, header=['Dose [Gy]', 'Relative Volume [%]'], index=False)
                DP_Histo_results.close()
                
                #Write DVH for BL plan
                BL_x = allVoxInt(BL_dose, struct) #This function calculates a 2D flat array of the dose for each voxel within the 3D structure
                BL_counts, BL_bin_edges = np.histogram(BL_x, bins=81, range=(0, 81), weights=None, density=False)
                BL_histcum = 100*(1 - np.cumsum(BL_counts)/len(BL_x)) #cumulative histogram values: y axis
                
                #Save DP_DVH to excel file
                BL_Histo_results = pd.ExcelWriter(BL_DVH_dir +'/' + file_write_name)
                BL_df_dict = {'Dose [Gy]': BL_bin_edges[:-1], 'Relative Volume [%]': BL_histcum}
                BL_Histo_df = pd.DataFrame(data=BL_df_dict)
                BL_Histo_df.to_excel(BL_Histo_results, header=['Dose [Gy]', 'Relative Volume [%]'], index=False)
                BL_Histo_results.close()
            
