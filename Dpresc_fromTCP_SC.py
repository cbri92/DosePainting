# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:31:51 2022

@author: Caterina Brighi

This script is used to take information about cellularity, alpha and beta and use them
to derive patient-specific dose prescriptions by means of a TCP model (poissonian, LQ).
"""


#%% Import functions 

import SimpleITK as sitk
import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from scipy import stats
import scipy.integrate as integrate
from ImageAnalysisFunctions import gunzip_shutil

def allVoxInt(image, roi):
    
    '''Returns a flatten array (in 2D) of the intensity values from all pixels in roi applied to an image.'''
    
    image = sitk.GetArrayFromImage(image)
    roi = sitk.GetArrayFromImage(roi)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    masked = image[~np.array(mask)]
    return masked

def generate_mask(image, roi):
    
    '''Returns the masked image of the roi applied to the image.
    Remember to save the masked image file after applying this function.'''
    
    masking = sitk.MaskImageFilter() 
    mask = masking.Execute(image, roi)
    return mask

def integrate_pdf(N):
    """This function calculates the integral of the kernel density estimation of the cellularity distribution.
        N: array of range of values of cellularity"""
    return integrate.quad(lambda N: kde(N), 0, N)[0]
    
def cdf_kde(X):
    """This function returns the cumulative density function of the cellularity distribution.
        X: array of range of values of cellularity"""
    if np.size(X) > 1:
        return np.array([integrate_pdf(x) for x in X])
    else:
        return integrate_pdf(X)

def transfer_func(X, Dmin, Dmax, N_dmax):
    """This function converts values of cellularity in an array X into dose values.
        Inputs:
            X: array of values of cellularity
            Dmin: minimum dose assigned to the CTV, i.e. 54 Gy
            Dmax: Dose corresponding to the 99th percentile of the patient-specific TCP
            N_dmax: mode value of cellularity.
           Outputs:
            Dose: array of dose
        """
    return (Dmax-Dmin)*cdf_kde(X)/cdf_kde(N_dmax) + Dmin

#%%Set alpha and beta particle-cell line specific values

# #Set alpha/beta photons
# alpha_beta_x = 2.4 #Gy for chordomas
# alpha_x = 0.1 #Gy-1
# beta_x = alpha_x/alpha_beta_x #Gy-2

#Values a and b carbon derived from plots in Fig 4 of paper Kato et al. International Journal of Radiation Oncology*Biology*Physics, 2011.
#Derivation in script Carbon_a_B_estimates.py

a = 0.950148844306803
a_std = 0.142982098449624
# b = 0.0677197997438524
# b_std = 0.0538530499689479
# a_b_carbon = a/b

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/nifti/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    MRI_dir = data_supradir+current+'/MRI/'

    studies = [ s.name for s in os.scandir(MRI_dir) if s.is_dir() ]
    studies = [x.replace('study_', '') for x in studies] #Remove word "study_" from list of studies
    studies2 = [int(x) for x in studies] #Convert list of strings in list of integers
    small_study = min(studies2) #Find the smallest study date
    study = 'study_'+str(small_study) #Define which one is the baseline (earliest) study

    if not os.path.exists(MRI_dir+study+'/micro_denoise/'):
        print('No images at baseline for '+current)
    else:
        
        subj_dir = MRI_dir+study+'/micro_denoise/' #Set path to subject directory where cell density map image and GTV image are stored
        if not os.path.exists(MRI_dir+study+'/Prescribed_dose'):#if it does not already exist, create a directory where the optimized dose will be saved
            os.mkdir(MRI_dir+study+'/Prescribed_dose')
        out_dir = MRI_dir+study+'/Prescribed_dose'
        
        if not os.path.exists(MRI_dir+study+'/micro_denoise/cellApp_CORRECTED.nii'):
            print('No cellApp image at baseline for '+current)
        else:
            #Unzip cellApp_CORRECTED.nii.gz file and delete original gzipped file
            for filename in glob.glob(MRI_dir+study+'/micro_denoise/cellApp_CORRECTED.nii.gz'):
                gunzip_shutil(filename, filename[:-3])
                
            #%%Generate CTV in DWI excluding bones
            
            CT_reg = sitk.ReadImage(MRI_dir+study+'/regImages/CTonDWI.nii')
            CTV_reg = sitk.ReadImage(MRI_dir+study+'/regImages/CTVonDWI.nii')
            bone_mask = CT_reg>100
            bone_mask = CTV_reg*bone_mask
            CTV_noBone = CTV_reg-bone_mask
            
            sitk.WriteImage(CTV_noBone, MRI_dir+study+'/micro_denoise/CTV_inDWI_noBone.nii')
            
            #%%
             
            print('Optimizing the dose prescription for '+current) 
                      
            cell_map = sitk.ReadImage(MRI_dir+study+'/micro_denoise/cellApp_CORRECTED.nii') #Read cellular density map. Units: um-3
            cell_map = cell_map*(10**8) #Convert cell map to Units: 10^4 cm-3
            
            cell_dir = MRI_dir+study+'/micro_denoise/' #Path to dir with noBone data
            CTV_full = sitk.ReadImage(cell_dir+'CTV_inDWI_noBone.nii') #Read CTV ROI without voxels including bones
            cell_map=generate_mask(cell_map, CTV_full) #Mask cellularity map only in CTV noBone
            sitk.WriteImage(cell_map, cell_dir+'cellApp_CORRECTED_noBone.nii') #Write cellularity map in valid voxels in noBone dir
            
            ctv_ok=sitk.ReadImage(MRI_dir+study+'/micro_denoise/ctv_voxel_ok_3D.nii.gz', sitk.sitkUInt8) #Read mask of CTV voxels that have valid values of ADC
            ctv_ok=(CTV_full+ctv_ok)>1 #Create an ROI of voxels in CTV that are valid and don't involve bone structures
            sitk.WriteImage(ctv_ok, cell_dir+'ctv_voxel_ok_3D_noBone.nii.gz')
            
            cell_type = 'noBone'
                               
            CTV = sitk.ReadImage(cell_dir+'ctv_voxel_ok_3D_noBone.nii.gz')
            CTV = sitk.Cast(CTV, sitk.sitkUInt8)
            CTV_outliers = CTV_full-CTV #Generate a mask of the outliers voxels by subtracting the voxels that have valid values of ADC from the CTV mask with noBone
            sitk.WriteImage(CTV_outliers, cell_dir+'CTV_outliers.nii')
            
            cell_map = generate_mask(cell_map, CTV) #Mask the cellularity map only in voxels of CTV that are valid for ADC and have no bones
            cell_map_nda = sitk.GetArrayFromImage(cell_map) #Convert cell map image into a numpy array  
            
            #Find indexes of voxels that have values of cellularity
            result = np.where(cell_map_nda>0)
            listOfCoordinates= list(zip(result[0], result[1], result[2])) #Save indexes in a list of coordinates 
        
            orig_shape = sitk.GetArrayFromImage(cell_map).shape #Save the shape of the original cell map image
            voxel_volume = (cell_map.GetSpacing()[0]*cell_map.GetSpacing()[1]*cell_map.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
            
            #Define input parameters for TCP modelling: number of cells per voxel (N0) and alpha and beta values per voxel
            CTV_final = CTV
            
            cell_dens = allVoxInt(cell_map, CTV_final) #Define array of values with cellular density>0
            N0 = cell_dens*voxel_volume #Multiply the values of cellular density per voxel by the volume of one voxel to obtain the number of cells per voxel
            
            alpha=np.full(len(N0), a, dtype=np.float32)
            # beta=np.full(len(N0), b, dtype=np.float32)
            
            # n=16 #Number of fractions in HD volume, while 9 fcn in LD volume
        
        #%%Find the patient specific Dmax corresponding to the 99th percentile of the TCP
            
            # for graphing below only
            x_vals = np.arange(0,100)
            #TCP Poisson model
            # vals=[np.prod( np.exp( - N0 * np.exp(-alpha*x-(beta*(x**2))/n) )) for x in x_vals]
            vals=[np.prod( np.exp( - N0 * np.exp(-alpha*x))) for x in x_vals]
            
            TCP = np.array(vals)
            VAL_max = np.isclose(TCP, 0.99, rtol=1e-02)
            MAX_DOSE = x_vals[VAL_max][0]
            print('MAX DOSE', MAX_DOSE)
            MIN_DOSE = x_vals[TCP>0.01][0]
            print('MIN DOSE', MIN_DOSE)
            
            # #Plot the TCP vs Dose
            # plt.plot(x_vals,TCP)
            # plt.plot([MAX_DOSE,MAX_DOSE],[0,1])
            # plt.plot([MIN_DOSE,MIN_DOSE],[0,1])
            # plt.legend(['TCP','Max dose value: '+str(MAX_DOSE), 'Min dose value: '+str(MIN_DOSE)], loc ="upper left")
            # plt.xlabel("Dose (Gy)")
            # plt.ylabel("TCP")
            # plt.savefig(out_dir+'/TCP_dose_plot_'+cell_type+'.jpg', bbox_inches='tight')
            # plt.show()
        
            # #Plot the probability density function vs dose
            # plt.plot(x_vals[:-1],np.diff(TCP))
            # plt.plot([MAX_DOSE,MAX_DOSE],[0,max(np.diff(TCP))])
            # plt.plot([MIN_DOSE,MIN_DOSE],[0,max(np.diff(TCP))])
            # plt.legend(['PDF','Max dose value: '+str(MAX_DOSE), 'Min dose value: '+str(MIN_DOSE)], loc ="upper left")
            # plt.xlabel("Dose (Gy)")
            # plt.ylabel("Probability density")
            # plt.savefig(out_dir+'/PDF_dose_plot_'+cell_type+'.jpg', bbox_inches='tight')
            # plt.show()
            
        # #%% Way of interpolating dose to obtain heterogeneous dose plan, such that to target higher dose to more radioresistant regions, whilst sparing dose to radiosensitive tissue (and surrounding healthy tissues)
        
        #     #%%Plot the probability density function vs n cells x vox
        
        #     plt.hist(N0, bins=200, density=True)
        #     plt.xlabel("N Cells per voxel")
        #     plt.ylabel("Probability density function")
        #     plt.savefig(out_dir+'/N_cells_histo_'+cell_type+'.jpg', bbox_inches='tight')
        #     plt.show()
        #     plt.close()
            
        #     #%%Plot the PDF and CDF
        #     fig, ax = plt.subplots()
        #     kde=stats.gaussian_kde(N0) #Estimate kernel density function
        #     X = np.linspace(0, max(N0)*1.1, 200)
        #     PDF = kde(X) #Calculate the PDF from the estimated kde, in the range of cellularity values available for this patient
        #     CDF = cdf_kde(X) #Calculate the cumulative density function from the estimated kde, in the range of cellularity values available for this patient
        #     ax2 = ax.twinx()
        #     ax2.fill_between(X, CDF, fc="#53e686", alpha=0.5)
        #     ax.fill_between(X, PDF, fc="#AAAAFF")
        #     ax.set_ylabel('Probability density function')
        #     ax2.set_ylabel('Cumulative density function')
        #     ax.set_xlabel('N of cells/voxel')
        #     fig.legend(['PDF', 'CDF'], loc=1)
        #     plt.show()
        #     fig.savefig(out_dir+'/PDF_CDF_N_'+cell_type+'_plot.jpg', bbox_inches='tight')
        #     plt.close()      
            
        #     #%%Map values of cellularity to values of dose
        #     Dmin, Dmax = 39.6, MAX_DOSE #Assign minimum and maximum values of dose to assign to voxels with lowest value of cellularity and mode (most frequent) values of cellularity, respectively. These for us are 39.6 Gy (corresponding to dose assigned to CTV_low dose) and MAX_DOSE, specific to the patient.
        #     pdf_vals = tuple(zip(kde(X), X))
        #     N_dmax = max(pdf_vals)[1]
        #     dose = transfer_func(N0, Dmin, Dmax, N_dmax)
            
        #     #Set values of dose > 80.96 Gy to 80.96 Gy as 110% of dose delivered to HD target volume in boost plans
        #     super_threshold_indices = dose > 80.96
        #     dose[super_threshold_indices] = 80.96
            
            
        #     #%%Plot the PDF and derived dose values as function of cellularity
        #     fig, ax = plt.subplots()
        #     plt.scatter(N0, dose, marker='x')
        #     ax2 = ax.twinx()
        #     ax.plot([N_dmax,N_dmax],[Dmin, max(dose)], 'k--')
        #     ax.plot([0,max(N0)],[Dmax, Dmax], 'k--')
        #     ax.annotate('Max Dose {0}'.format(Dmax), (max(N0)*0.7, Dmax*1.02))
        #     ax.annotate('mode of N density {0}'.format(int(N_dmax)), (N_dmax*1.1, Dmin*1.02))
        #     ax2.fill_between(X, PDF, fc="#AAAAFF", alpha=0.5)
        #     ax.set_ylabel('Dose [Gy]')
        #     ax2.set_ylabel('Probability density function')
        #     ax.set_xlabel('N of cells/voxel')
        #     fig.savefig(out_dir+'/PDF_DOSE_N_'+cell_type+'_plot.jpg', bbox_inches='tight')
        #     plt.show()
        #     plt.close()
            
        #     #%% Compare new TCP with TCP that would be obtained by deliverying an homogeneous dose boost to the GTV to maximum dose estimated in the optimization
        #     # TCP_heter =np.prod( np.exp( - N0 * np.exp(-alpha*dose-(beta*(dose**2))/n))) #TCP obtained with optimised heterogeneous dose distribution
        #     # TCP_boost = np.prod( np.exp( - N0 * np.exp(-alpha*Dmax-(beta*(Dmax**2))/n))) # TCP with optimised homogeneous max dose boost to entire GTV
        #     TCP_heter =np.prod( np.exp( - N0 * np.exp(-alpha*dose))) #TCP obtained with optimised heterogeneous dose distribution
        #     TCP_boost = np.prod( np.exp( - N0 * np.exp(-alpha*Dmax)))# TCP with optimised homogeneous max dose boost to entire GTV
        #     print('TCP optimised heterogeneous dose plan: ',TCP_heter)
        #     print('TCP optimised homogeneous max dose boost plan: ',TCP_boost)
            
        #     plt.hist(dose, bins=20)
        #     plt.xlabel("Dose (Gy)")
        #     plt.ylabel("Counts")
        #     plt.savefig(out_dir+'/Dose_histo_'+cell_type+'.jpg', bbox_inches='tight')
        #     plt.show()
        #     plt.close()
            
        # #%% Saving the optimised dose values into a prescription dose image
            
        #     dose_nda = np.ndarray(shape=orig_shape) #Generate a 3D array of the same shape of the cell density map
            
        #     for d,i in zip(list(dose),listOfCoordinates): #Assign the calculated optimised dose to each voxel, based on the index of the corresponding cellularity voxel
        #         dose_nda[i]=d
        
        #     dose_img = sitk.GetImageFromArray(dose_nda) #Convert the numpy array into an image 
        #     dose_img.CopyInformation(cell_map) #Copy the headers of the cellularity image onto the headers of the optimised dose image
            
        #     sitk.WriteImage(dose_img, out_dir+'/Dose_optimised_CTV'+cell_type+'.nii') #Save the prescribed optimised dose image as a nifti file 
            
        # #%%Set values within CTV with cellularity untrusted and bone pixels to receive Dmin = 54 Gy
            
        #     # CTV_complete = sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/CTV_inDWI_ITK.nii')
        #     # CTV_remaining = (CTV_complete>0) & (dose_img==0)
        #     # dose_img_final = set_mask_value(dose_img, CTV_remaining, Dmin)
        #     # sitk.WriteImage(dose_img_final, out_dir+'/Dose_optimised_CTV'+cell_type+'_final54.nii')
            