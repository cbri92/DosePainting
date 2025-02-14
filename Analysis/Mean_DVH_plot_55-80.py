#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:03:28 2024

@author: Caterina Brighi
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%% Set Working directory
        
data_supradir = '/Users/cbri3325/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)(2)/Caterina Brighi/Data/SBC_Tutti/' #Set working directory
DVH_dir = data_supradir+'Volumes_DVH/'

#%%Read GTV dataframe from excel spreadsheets and plotting DVH plot

print('Reading GTV dataframe from excel spreadsheet and plotting DVH')

GTV_file = DVH_dir+'GTV_DVH.xlsx'

GTV_DP_df = pd.read_excel(
    GTV_file,
    sheet_name='Dose painting',  # specific sheet
    usecols=['Dose [Gy]', 'Mean Volume [%]','25th percentile Volume [%]','75th percentile Volume [%]'],  # specific columns
    skiprows=range(1,56),  # skip the first 3 rows
    header=0  # use the first row as column names
)

GTV_BL_df = pd.read_excel(
    GTV_file,
    sheet_name='Baseline',  # specific sheet
    usecols=['Dose [Gy]', 'Mean Volume [%]','25th percentile Volume [%]','75th percentile Volume [%]'],  # specific columns
    skiprows=range(1,56),  # skip the first 3 rows
    header=0  # use the first row as column names
)

#Generate DVH plots
colors = ["#632de9", "#8e82fe", "#13eac9", "#7bf2da", "#f97306", "#ffb07c", "#15b01a", "#6fc276"] 
plt.plot(GTV_DP_df['Dose [Gy]'], GTV_DP_df['Mean Volume [%]'], label='Dose painting', color=colors[6], linewidth=3)
plt.fill_between(GTV_DP_df['Dose [Gy]'], GTV_DP_df['25th percentile Volume [%]'], GTV_DP_df['75th percentile Volume [%]'],alpha=0.4, edgecolor=colors[6], facecolor=colors[7])
plt.plot(GTV_BL_df['Dose [Gy]'], GTV_BL_df['Mean Volume [%]'], label='Uniform', color=colors[4], linewidth=3)
plt.fill_between(GTV_BL_df['Dose [Gy]'], GTV_BL_df['25th percentile Volume [%]'], GTV_BL_df['75th percentile Volume [%]'],alpha=0.4, edgecolor=colors[4], facecolor=colors[5])
plt.legend(loc='best', fontsize=15)
plt.xlabel('Dose [Gy]',fontsize=15)
plt.ylabel('Volume [%]',fontsize=15)
plt.xlim(55,80)
plt.ylim(0,100)
plt.xticks([55,60,65,70,75,80])
plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.title('GTV',fontsize=15)     
plt.grid(color = 'grey', linestyle = '--', linewidth = 0.5, alpha=0.5)
plt.tight_layout()
plt.savefig(DVH_dir +'GTV_DVH_55-80.png')
plt.show()
plt.close()


#%%Read CTV dataframe from excel spreadsheets and plotting DVH plot

print('Reading CTV dataframe from excel spreadsheet and plotting DVH')

CTV_file = DVH_dir+'CTV_DVH.xlsx'

CTV_DP_df = pd.read_excel(
    CTV_file,
    sheet_name='Dose painting',  # specific sheet
    usecols=['Dose [Gy]', 'Mean Volume [%]','25th percentile Volume [%]','75th percentile Volume [%]'],  # specific columns
    skiprows=range(1,56),  # skip the first 3 rows
    header=0  # use the first row as column names
)

CTV_BL_df = pd.read_excel(
    CTV_file,
    sheet_name='Baseline',  # specific sheet
    usecols=['Dose [Gy]', 'Mean Volume [%]','25th percentile Volume [%]','75th percentile Volume [%]'],  # specific columns
    skiprows=range(1,56),  # skip the first 3 rows
    header=0  # use the first row as column names
)

#Generate DVH plots
colors = ["#632de9", "#8e82fe", "#13eac9", "#7bf2da", "#f97306", "#ffb07c", "#15b01a", "#6fc276"] 
plt.plot(CTV_DP_df['Dose [Gy]'], CTV_DP_df['Mean Volume [%]'], label='Dose painting', color=colors[6], linewidth=3)
plt.fill_between(CTV_DP_df['Dose [Gy]'], CTV_DP_df['25th percentile Volume [%]'], CTV_DP_df['75th percentile Volume [%]'],alpha=0.4, edgecolor=colors[6], facecolor=colors[7])
plt.plot(CTV_BL_df['Dose [Gy]'], CTV_BL_df['Mean Volume [%]'], label='Uniform', color=colors[4], linewidth=3)
plt.fill_between(CTV_BL_df['Dose [Gy]'], CTV_BL_df['25th percentile Volume [%]'], CTV_BL_df['75th percentile Volume [%]'],alpha=0.4, edgecolor=colors[4], facecolor=colors[5])
plt.legend(loc='best', fontsize=15)
plt.xlabel('Dose [Gy]',fontsize=15)
plt.ylabel('Volume [%]',fontsize=15)
plt.xlim(55,80)
plt.ylim(0,100)
plt.xticks([55,60,65,70,75,80])
plt.yticks([0,10,20,30,40,50,60,70,80,90,100])
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.title('CTV',fontsize=15)     
plt.grid(color = 'grey', linestyle = '--', linewidth = 0.5, alpha=0.5)
plt.tight_layout()
plt.savefig(DVH_dir +'CTV_DVH_55-80.png')
plt.show()
plt.close()