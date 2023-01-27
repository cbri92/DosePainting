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
        
# data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO' #Set working directory
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/test/nifti2/'

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current #Define path to patient directory
    
    CT = sitk.ReadImage(subj_dir+'/ct.nii', sitk.sitkFloat32) #Read CT image
    #Aggiungi trova GTV, CTV, PTV  
    
    MRI_dir = subj_dir+'/MRI/'
    studies_path = [ f.path for f in os.scandir(MRI_dir) if f.is_dir() ] #Create a list of the paths to the MRI studies directories
    studies_name = [ f.name for f in os.scandir(MRI_dir) if f.is_dir() ] #Create a list of studies names
    
    for study in studies_name:
        
        study_dir = MRI_dir+study
        
        if not os.listdir(study_dir) :
            print("Study "+study+" directory is empty")
        else:           
            #Create a directory to store the registration matrices
            if not os.path.exists(study+'/regFiles'):#if it does not already exist, create a directory where the registration files will be saved
                os.mkdir(study_dir+'/regFiles')
            reg_dir = study_dir+'/regFiles'
            
            #Create a directory to store the registered images
            if not os.path.exists(study+'/regImages'):#if it does not already exist, create a directory where the registered images will be saved
                os.mkdir(study_dir+'/regImages')
            regIm_dir = study_dir+'/regImages'
            
            #Create a directory to store the MRI info files
            if not os.path.exists(study_dir+'/info/'):#if it does not already exist, create a directory where the info files will be saved
                os.mkdir(study_dir+'/info/')
            info_dir = study_dir+'/info/'
            
            #Move all the .bval, .bvec and .json files into the info folder
            for filename in glob.glob(study_dir +'/' +'*.bval'):
                shutil.move(filename, info_dir)
                
            for filename in glob.glob(study_dir +'/' +'*.bvec'):
                shutil.move(filename, info_dir)
                
            for filename in glob.glob(study_dir +'/' +'*.json'):
                shutil.move(filename, info_dir)
                
            #%%Registration of CT image to T2 MRI image
            
            #Read T2 transversal image       
            T2_path = glob.glob(study_dir+'/'+'*t2_tse_tra*')[0]
            T2 = sitk.ReadImage(T2_path, sitk.sitkFloat32)
            
            #Apply rigid body registration to register CT image to T2 image
            CTinT2, CTinT2_tfm = RegisterResample_image(T2, CT, 'different', 'Multires')
            sitk.WriteImage(CTinT2, regIm_dir+'/CTonT2.nii') #Save registered image in regImages directory
            sitk.WriteTransform(CTinT2_tfm, reg_dir+'/xf_CTonT2.txt') #Save registration matrix in regFiles directory
            
            #%%Registration of T2 MRI image to DWI MRI image
            
            #Read DWI image
            for filename in glob.glob(study_dir+'/'+'*ep2d_diff*'):
                if (('nii' in filename) and ('adc' not in filename)):
                    DWI_path = filename
            
            DWI = sitk.ReadImage(DWI_path, sitk.sitkFloat32)
            
            #Apply rigid body registration to register T2 image to DWI image
            T2inDWI, T2inDWI_tfm = RegisterResample_image(DWI, T2, 'same', 'Multires')
            sitk.WriteImage(T2inDWI, regIm_dir+'/T2onDWI.nii') #Save registered image in regImages directory
            sitk.WriteTransform(T2inDWI_tfm, reg_dir+'/xf_T2onDWI.txt') #Save registration matrix in regFiles directory
