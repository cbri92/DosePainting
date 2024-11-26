#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:03:28 2024

@author: Caterina Brighi

This script generates mean DVH plots in the RT targets and in the organs at risk for the entire patients cohort for baseline and dose painting plans. The plots also contain the 25th and 75th percentile bands of the Volume %.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%% Set Working directory
        
data_supradir = 'path/to/patients/data/supra/directory/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

if not os.path.exists(data_supradir+'Volumes_DVH/'):
    os.mkdir(data_supradir+'Volumes_DVH/')
out_dir = data_supradir+'Volumes_DVH/'
 
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

structures = ['GTV','CTV','Brainstem','Chiasm','Optic nerve L','Optic nerve R','Coclea L','Coclea R','Temporal lobe L','Temporal lobe R','Carotid L','Carotid R']
PLANS = ['DP', 'BL']

for struc in structures:
    globals()[struc+'_DVH'] = pd.ExcelWriter(out_dir+struc+'_DVH.xlsx')
    globals()[struc+'_DP_df'] = pd.DataFrame()
    globals()[struc+'_BL_df'] = pd.DataFrame()
    
#%%Create a for loop to perform image analysis on each subject sequentially

for current,i in zip(subjs_name,range(len(subjs_name))):

    subj_dir = data_supradir+current #Set path to subject directory
    
    DP_DVH_dir = subj_dir + '/DVH_files2/Dose_painting_4B/'
    BL_DVH_dir = subj_dir + '/DVH_files2/Baseline_4B/'
    
    DVH_dirs = [DP_DVH_dir, BL_DVH_dir]   
    #%%Create lists of structures file paths and structures names
    DP_files_paths = []
    DP_structures_names = []
    
    BL_files_paths = []
    BL_structures_names = []
    
    for plan_dir,plan in zip(DVH_dirs, PLANS):
        
        files_paths = []
        structures_names = []
        
        print('Finding structures DVH for '+plan+' plan for '+current)
        
        for f in os.scandir(plan_dir):
            if f.name == "GTV.xlsx":
                print(f.name+' found')
                GTV_f = f.path
                GTV_name = "GTV"
                files_paths.append(GTV_f)
                structures_names.append(GTV_name)
                
            elif f.name == "CTV_AIRC24946.xlsx":
                print(f.name+' found')
                CTV_f = f.path
                CTV_name = "CTV"
                files_paths.append(CTV_f)
                structures_names.append(CTV_name)
                
            elif f.name in tronco:
                print(f.name+' found')
                brainstem_f = f.path
                brainstem_name = "Brainstem"
                files_paths.append(brainstem_f)
                structures_names.append(brainstem_name)
                
            elif f.name in chiasma:
                print(f.name+' found')
                chiasm_f = f.path
                chiasm_name = "Chiasm"
                files_paths.append(chiasm_f)
                structures_names.append(chiasm_name)
                
            elif f.name in nervo_ottico_sx:
                print(f.name+' found')
                opt_nerv_l_f = f.path
                opt_nerv_l_name = "Optic nerve L"
                files_paths.append(opt_nerv_l_f)
                structures_names.append(opt_nerv_l_name)
                
            elif f.name in nervo_ottico_dx:
                print(f.name+' found')
                opt_nerv_r_f = f.path
                opt_nerv_r_name = "Optic nerve R"
                files_paths.append(opt_nerv_r_f)
                structures_names.append(opt_nerv_r_name)
                
            elif f.name in coclea_sx:
                print(f.name+' found')
                coclea_l_f = f.path
                coclea_l_name = "Coclea L"
                files_paths.append(coclea_l_f)
                structures_names.append(coclea_l_name)
                
            elif f.name in coclea_dx:
                print(f.name+' found')
                coclea_r_f = f.path
                coclea_r_name = "Coclea R"
                files_paths.append(coclea_r_f)
                structures_names.append(coclea_r_name)
                
            elif f.name in temp_lobe_sx:
                print(f.name+' found')
                temp_lobe_l_f = f.path
                temp_lobe_l_name = "Temporal lobe L"
                files_paths.append(temp_lobe_l_f)
                structures_names.append(temp_lobe_l_name)
                
            elif f.name in temp_lobe_dx:
                print(f.name+' found')
                temp_lobe_r_f = f.path
                temp_lobe_r_name = "Temporal lobe R"
                files_paths.append(temp_lobe_r_f)
                structures_names.append(temp_lobe_r_name)
                
            elif f.name in carotide_sx:
                print(f.name+' found')
                carotid_l_f = f.path
                carotid_l_name = "Carotid L"
                files_paths.append(carotid_l_f)
                structures_names.append(carotid_l_name)
                
            elif f.name in carotide_dx:
                print(f.name+' found')
                carotid_r_f = f.path
                carotid_r_name = "Carotid R"
                files_paths.append(carotid_r_f)
                structures_names.append(carotid_r_name)
                
        if plan_dir is DP_DVH_dir:
        
            DP_files_paths = files_paths
            DP_structures_names = structures_names
            
        elif plan_dir is BL_DVH_dir:
            
            BL_files_paths = files_paths
            BL_structures_names = structures_names
     
    #%%Read Dose and Volume % data from each excel file and generate Dose %. Append to new dataframe Dose% and Volume% for each structure.
    paths_dirs = [DP_files_paths, BL_files_paths]
    names_struc_dirs = [DP_structures_names, BL_structures_names]
    
    print('Appending structures DVH to study dataframe for '+current)
    
    for pdir, nsdir, plan in zip(paths_dirs, names_struc_dirs, PLANS):
        
        for fpath, name in zip(pdir,nsdir):
            
            df=pd.read_excel(fpath, sheet_name='Sheet1', header=0)
            Volume_perc = np.array(df['Relative Volume [%]'])
            Dose = np.array(df['Dose [Gy]'])
            # if name in ['GTV','CTV']:
            #     Max_dose = 74
            # elif name == 'Brainstem':
            #     Max_dose = 61
            # elif name == 'Chiasm':
            #     Max_dose = 54
            # elif name in ['Optic nerve L', 'Optic nerve R']:
            #     Max_dose = 50
            # elif name in ['Coclea L', 'Coclea R']:
            #     Max_dose = 45
            # elif name in ['Temporal lobe L', 'Temporal lobe R']:
            #     Max_dose = 71
            # elif name in ['Carotid L', 'TCarotid R']:
            #     Max_dose = 75.5
            # # Max_dose = df[0].max()
            # Dose_perc = np.array((Dose/Max_dose)*100)
            
            if int(i)+1 == 1:
                globals()[name+'_'+plan+'_df']['Dose [Gy]'] = Dose
            globals()[name+'_'+plan+'_df']['P'+str(int(i)+1)+' Volume [%]'] = Volume_perc
         
            # str_data={'Dose [%]':Dose_perc, 'P'+str(int(i)+1)+' Volume [%]':Volume_perc}
            # globals()[name+'_'+plan+'_df'] = pd.concat([globals()[name+'_'+plan+'_df'], pd.DataFrame(str_data)], ignore_index=True)
            # globals()[name+'_'+plan+'_df'] = pd.concat([globals()[name+'_'+plan+'_df'], pd.from_dict(str_data, orient='columns'), ignore_index=True)


#%%Calculate average and 25th and 75th percentiles of Volumes % data

for struc in structures:
    
    for plan in PLANS:
        
        # String to search for
        search_string = ' Volume [%]'
        
        # Get columns whose names contain the search string
        columns_with_string = [col for col in globals()[struc+'_'+plan+'_df'].columns if search_string in col]
        
        print(current+' has '+str(columns_with_string)+' for '+struc+' '+plan+' plan')
        
        globals()['mean_'+plan] = globals()[struc+'_'+plan+'_df'][columns_with_string].mean(axis=1)
        globals()['p25_'+plan] =globals()[struc+'_'+plan+'_df'][columns_with_string].quantile(q=0.25,axis=1)
        globals()['p75_'+plan] =globals()[struc+'_'+plan+'_df'][columns_with_string].quantile(q=0.75,axis=1)
        
        globals()[struc+'_'+plan+'_df']['Mean Volume [%]'] = globals()['mean_'+plan]
        globals()[struc+'_'+plan+'_df']['25th percentile Volume [%]'] = globals()['p25_'+plan]
        globals()[struc+'_'+plan+'_df']['75th percentile Volume [%]'] = globals()['p75_'+plan]

#%%Save dataframes to excel spreadsheets and generate DVH plots

print('Saving study dataframes to excel spreadsheets and plot structures DVH')

for STRUCT in structures:
    
    globals()[STRUCT+'_DP_df'].to_excel(globals()[STRUCT+'_DVH'], sheet_name='Dose painting', index=False)
    globals()[STRUCT+'_BL_df'].to_excel(globals()[STRUCT+'_DVH'], sheet_name='Baseline', index=False)
    globals()[STRUCT+'_DVH'].close()   
    
    #Set percentile values to be equal to mean values where 25th percentile is >= mean, and where 75th percentile is < mean
    for plan in PLANS:
        globals()[STRUCT+'_'+plan+'_df']['25th percentile Volume [%]'] = globals()[STRUCT+'_'+plan+'_df']['25th percentile Volume [%]'].where(globals()[STRUCT+'_'+plan+'_df']['25th percentile Volume [%]'] <= globals()[STRUCT+'_'+plan+'_df']['Mean Volume [%]'], globals()[STRUCT+'_'+plan+'_df']['Mean Volume [%]'])
        globals()[STRUCT+'_'+plan+'_df']['75th percentile Volume [%]'] = globals()[STRUCT+'_'+plan+'_df']['75th percentile Volume [%]'].where(globals()[STRUCT+'_'+plan+'_df']['75th percentile Volume [%]'] > globals()[STRUCT+'_'+plan+'_df']['Mean Volume [%]'], globals()[STRUCT+'_'+plan+'_df']['Mean Volume [%]'])
    
    #Generate DVH plots
    colors = ["#632de9", "#8e82fe", "#13eac9", "#7bf2da", "#f97306", "#ffb07c", "#15b01a", "#6fc276"] 
    plt.plot(globals()[STRUCT+'_DP_df']['Dose [Gy]'], globals()[STRUCT+'_DP_df']['Mean Volume [%]'], label='Dose painting', color=colors[6], linewidth=3)
    plt.fill_between(globals()[STRUCT+'_DP_df']['Dose [Gy]'], globals()[STRUCT+'_DP_df']['25th percentile Volume [%]'], globals()[STRUCT+'_DP_df']['75th percentile Volume [%]'],alpha=0.4, edgecolor=colors[6], facecolor=colors[7])
    plt.plot(globals()[STRUCT+'_BL_df']['Dose [Gy]'], globals()[STRUCT+'_BL_df']['Mean Volume [%]'], label='Baseline', color=colors[4], linewidth=3)
    plt.fill_between(globals()[STRUCT+'_BL_df']['Dose [Gy]'], globals()[STRUCT+'_BL_df']['25th percentile Volume [%]'], globals()[STRUCT+'_BL_df']['75th percentile Volume [%]'],alpha=0.4, edgecolor=colors[4], facecolor=colors[5])
    plt.legend(loc='best', fontsize=15)
    plt.xlabel('Dose [Gy]',fontsize=15)
    plt.ylabel('Volume [%]',fontsize=15)
    plt.xlim(0,80)
    plt.ylim(0,100)
    plt.xticks([0,10,20,30,40,50,60,70,80])
    plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    plt.title(STRUCT,fontsize=15)     
    plt.grid(color = 'grey', linestyle = '--', linewidth = 0.5, alpha=0.5)
    plt.tight_layout()
    plt.savefig(out_dir + STRUCT+'_DVH.png')
    plt.show()
    plt.close()
