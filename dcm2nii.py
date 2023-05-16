# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 10:33:42 2022

@author: Caterina Brighi
"""

import os
import glob
import shutil
import SimpleITK as sitk
import json
import sys
import dicom2nifti
sys.path.append("C:/Users/cbri3325/Anaconda3/lib/site-packages/platipy/__init__.py") # Path containing PlatiPy library
from platipy.dicom.io.rtstruct_to_nifti import convert_rtstruct
from platipy.dicom.io.rtdose_to_nifti import convert_rtdose
from platipy.dicom.io.rtstruct_to_nifti import read_dicom_struct_file


def DICOMseries_toNII(dcm_dir_path, nii_filepath, headers_filepath):
    '''This function converts DICOM series to nii images and saves also a json file with the dicom headers.'''
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(dcm_dir_path)
    reader.SetFileNames(dicom_names)
    image = reader.Execute()
    sitk.WriteImage(image, nii_filepath+'.nii')
    
    reader1 = sitk.ImageFileReader()
    reader1.SetFileName(dicom_names[0])
    reader1.LoadPrivateTagsOn()
    reader1.ReadImageInformation()    
    headers ={}
    for k in reader1.GetMetaDataKeys():
        v = reader1.GetMetaData(k)
        headers[k]=v
        # print(f'({k}) = "{v}"')
    
    with open(headers_filepath+'.json', 'w') as fp:
        json.dump(headers, fp)

    
nii_dir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/nifti/'
MRI_dir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/dicom/MRI/'
RT_dir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/dicom/RTDATA/'

#%%Convert MRI data

print('Converting MRI data from dicom to nifti...')
subjs_path = [ f.path for f in os.scandir(MRI_dir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(MRI_dir) if f.is_dir() ] #Create a list of subjects names

for current in subjs_name:
    
    print('for patient '+current)

    nifti_path_pt = os.path.join(nii_dir, current)
    if not os.path.exists(nifti_path_pt):
        os.mkdir(nifti_path_pt)
    nii_MRI_dir = os.path.join(nifti_path_pt, 'MRI')
    if not os.path.exists(nii_MRI_dir):
        os.mkdir(nii_MRI_dir)
        
    subj_dir = MRI_dir+current #Set path to subject directory    
    
    studies_path = [ f.path for f in os.scandir(subj_dir) if f.is_dir() ] #Create a list of the paths to the subjects directories
    studies_name = [ f.name for f in os.scandir(subj_dir) if f.is_dir() ]
    
    for study in studies_name:
        MRI_study_dir = os.path.join(nii_MRI_dir, study)
        if not os.path.exists(MRI_study_dir):
            os.mkdir(MRI_study_dir)
        
        study_dir = os.path.join(subj_dir, study+'/MR/')
        scans_path = [ f.path for f in os.scandir(study_dir) if f.is_dir() ]
        scans_name = [ f.name for f in os.scandir(study_dir) if f.is_dir() ]
        
        for path, name in zip(scans_path, scans_name):
            # DICOMseries_toNII(path, os.path.join(MRI_study_dir, name), os.path.join(MRI_study_dir, name))

            # scan_nii_dir = os.path.join(MRI_study_dir, name)
            # if not os.path.exists(scan_nii_dir):
            #     os.mkdir(scan_nii_dir)
            dicom2nifti.convert_directory(path,MRI_study_dir)
        
#%% Convert RT data    

print('Converting RT data from dicom to nifti...')
subjs_path = [ f.path for f in os.scandir(RT_dir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(RT_dir) if f.is_dir() ] #Create a list of subjects names

for current in subjs_name:
    
    print('for patient '+current)

    nifti_path_pt = os.path.join(nii_dir, current)
    if not os.path.exists(nifti_path_pt):
        os.mkdir(nifti_path_pt)
    
    nii_RTDOSE_dir = os.path.join(nifti_path_pt, 'RTDOSE')
    if not os.path.exists(nii_RTDOSE_dir):
        os.mkdir(nii_RTDOSE_dir)
        
    nii_RTSTRUCT_dir = os.path.join(nifti_path_pt, 'RTSTRUCT')
    if not os.path.exists(nii_RTSTRUCT_dir):
        os.mkdir(nii_RTSTRUCT_dir)
    
    subj_dir = RT_dir+current #Set path to subject directory 
    
    CT_dir = subj_dir+'/CT/'
    DICOMseries_toNII(CT_dir, nifti_path_pt+'/ct', nifti_path_pt+'/ct')
    # dicom2nifti.convert_dicom.dicom_series_to_nifti(CT_dir,nifti_path_pt+'/ct.nii') 
                  
    RTSTRUCT_dir = subj_dir+'/RTSTRUCT/'
    PHYS_dir = subj_dir+'/PHYS/'
    LET_dir = subj_dir+'/LET/'
    RBE_dir = subj_dir+'/LEM/'
    
    #Check and correct any label in RTSTRUCT that has '/' in the naming
    rtstr = read_dicom_struct_file(RTSTRUCT_dir+'RTSTRUCT.dcm')
    for i in rtstr.StructureSetROISequence:
        ch = '/'
        if ch in i.ROIName:
            x = i.ROIName.replace(ch,'_')
            i.ROIName = x
    rtstr.save_as(RTSTRUCT_dir+'RTSTRUCT.dcm')
    
    #Convert RTSTRUCT
    convert_rtstruct(CT_dir, RTSTRUCT_dir+'RTSTRUCT.dcm', prefix='', output_dir=nii_RTSTRUCT_dir)
    
    #Convert Dose files
    convert_rtdose(PHYS_dir+'FROG_PHYS_1.dcm', False, nii_RTDOSE_dir+'/PHYS_v0.nii')
    convert_rtdose(PHYS_dir+'FROG_PHYS_2.dcm', False, nii_RTDOSE_dir+'/PHYS_v1.nii')
    convert_rtdose(LET_dir+'FROG_LETd_1.dcm', False, nii_RTDOSE_dir+'/LET1.nii')
    convert_rtdose(LET_dir+'FROG_LETd_2.dcm', False, nii_RTDOSE_dir+'/LET2.nii')
    convert_rtdose(RBE_dir+'FROG_LEM_1.dcm', False, nii_RTDOSE_dir+'/RBE_v0.nii')
    convert_rtdose(RBE_dir+'FROG_LEM_2.dcm', False, nii_RTDOSE_dir+'/RBE_v1.nii')
        
    
    