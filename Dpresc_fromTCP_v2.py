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
from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import *
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
        """This function converts values of cellularity in an aaray X into dose values.
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

#Set alpha/beta photons
alpha_beta_x = 2.4 #Gy for chordomas
alpha_x = 0.1 #Gy-1
beta_x = alpha_x/alpha_beta_x #Gy-2

#Values RBEmax and RBEmin for proton on chordoma cells taken from paper Paganetti, International Journal of Radiation Oncology*Biology*Physics, 2022, 112(1), 222-236.

RBEmax = 1.59
RBEmin = 1.18

a = (RBEmax*alpha_x)
b = (beta_x*(RBEmin**2))
#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('AIRC24946_R052')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/MRI/baseline/micro/' #Set path to subject directory where cell density map image and GTV image are stored
    if not os.path.exists(data_supradir+current+'/MRI/baseline/Prescribed_dose'):#if it does not already exist, create a directory where the optimized dose will be saved
        os.mkdir(data_supradir+current+'/MRI/baseline/Prescribed_dose')
    out_dir = data_supradir+current+'/MRI/baseline/Prescribed_dose'
    
    #Unzip cellApp_CORRECTED.nii.gz file and delete original gzipped file
    for filename in glob.glob(data_supradir+current+'/MRI/baseline/orig/' +'cellApp_CORRECTED.nii.gz'):
        gunzip_shutil(filename, filename[:-3])
     
    print('Optimizing the dose prescription for '+current) 
              
    cell_map = sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/cellApp_CORRECTED.nii') #Read cellular density map. Units: um-3
    cell_map = cell_map*(10**8) #Convert cell map to Units: 10^4 cm-3
    
    cell_dirs = [data_supradir+current+'/MRI/baseline/orig/',data_supradir+current+'/MRI/baseline/noBone/']
    cell_types = ['orig', 'noBone']
    
    for cell_dir, cell_type in zip(cell_dirs,cell_types):

        if cell_dir == data_supradir+current+'/MRI/baseline/orig/':
            CTV_full = sitk.ReadImage(cell_dir+'CTV_inDWI_ITK.nii') #Read CTV
            
        elif cell_dir == data_supradir+current+'/MRI/baseline/noBone/':
            CTV_full = sitk.ReadImage(cell_dir+'CTV_inDWI_noBone.nii') #Read CTV
            
            cell_map=sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/cellApp_CORRECTED.nii')
            cell_map=generate_mask(cell_map, CTV_full)
            cell_map = cell_map*(10**8)
            sitk.WriteImage(cell_map, cell_dir+'cellApp_CORRECTED.nii')
            
            ctv_ok=sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/ctv_voxel_ok_3D.nii.gz', sitk.sitkUInt8)
            ctv_ok=(CTV_full+ctv_ok)>1
            sitk.WriteImage(ctv_ok, cell_dir+'ctv_voxel_ok_3D.nii.gz')
                       
        CTV = sitk.ReadImage(cell_dir+'ctv_voxel_ok_3D.nii.gz')
        CTV = sitk.Cast(CTV, sitk.sitkUInt8)
        CTV_outliers = CTV_full-CTV
        sitk.WriteImage(CTV_outliers, cell_dir+'CTV_outliers.nii')
        
        cell_map = generate_mask(cell_map, CTV)
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
        beta=np.full(len(N0), b, dtype=np.float32)
        
        n=37 #Number of fractions
    
    #%%Find the patient specific Dmax corresponding to the 99th percentile of the TCP
        
        # for graphing below only
        x_vals = np.arange(0,100)
        #TCP Poisson model
        vals=[np.prod( np.exp( - N0 * np.exp(-alpha*x-(beta*(x**2))/n) )) for x in x_vals]
        
        TCP = np.array(vals)
        VAL_max = np.isclose(TCP, 0.99, rtol=1e-02)
        MAX_DOSE = x_vals[VAL_max][0]
        print('MAX DOSE', MAX_DOSE)
        MIN_DOSE = x_vals[TCP>0.01][0]
        print('MIN DOSE', MIN_DOSE)
        
        #Plot the TCP vs Dose
        plt.plot(x_vals,TCP)
        plt.plot([MAX_DOSE,MAX_DOSE],[0,1])
        plt.plot([MIN_DOSE,MIN_DOSE],[0,1])
        plt.legend(['TCP','Max dose value: '+str(MAX_DOSE), 'Min dose value: '+str(MIN_DOSE)], loc ="upper left")
        plt.xlabel("Dose (Gy)")
        plt.ylabel("TCP")
        plt.savefig(out_dir+'/TCP_dose_plot_'+cell_type+'.jpg')
        plt.show()
    
        #Plot the probability density function vs dose
        plt.plot(x_vals[:-1],np.diff(TCP))
        plt.plot([MAX_DOSE,MAX_DOSE],[0,max(np.diff(TCP))])
        plt.plot([MIN_DOSE,MIN_DOSE],[0,max(np.diff(TCP))])
        plt.legend(['PDF','Max dose value: '+str(MAX_DOSE), 'Min dose value: '+str(MIN_DOSE)], loc ="upper left")
        plt.xlabel("Dose (Gy)")
        plt.ylabel("Probability density")
        plt.savefig(out_dir+'/PDF_dose_plot_'+cell_type+'.jpg')
        plt.show()
        
    #%% Way of interpolating dose to obtain heterogeneous dose plan, such that to target higher dose to more radioresistant regions, whilst sparing dose to radiosensitive tissue (and surrounding healthy tissues)
    
        #%%Plot the probability density function vs n cells x vox

        plt.hist(N0, bins=200, density=True)
        plt.xlabel("N Cells per voxel")
        plt.ylabel("Probability density function")
        plt.savefig(out_dir+'/N_cells_histo_'+cell_type+'.jpg')
        plt.close()
        
        #%%Plot the PDF and CDF
        fig, ax = plt.subplots()
        kde=stats.gaussian_kde(N0) #Estimate kernel density function
        X = np.linspace(0, max(N0)*1.1, 200)
        PDF = kde(X) #Calculate the PDF from the estimated kde, in the rance of cellularity values available for this patient
        CDF = cdf_kde(X) #Calculate the cumulative density function from the estimated kde, in the rance of cellularity values available for this patient
        ax2 = ax.twinx()
        ax2.fill_between(X, CDF, fc="#53e686", alpha=0.5)
        ax.fill_between(X, PDF, fc="#AAAAFF")
        ax.set_ylabel('Probability density function')
        ax2.set_ylabel('Cumulative density function')
        ax.set_xlabel('N of cells/voxel')
        fig.legend(['PDF', 'CDF'], loc=1)
        plt.show()
        plt.savefig(out_dir+'/PDF_CDF_N_'+cell_type+'_plot.jpg')
        plt.close()
        
        
        #%%Map values of cellularity to values of dose
        Dmin, Dmax = 54, MAX_DOSE #Assign minimum and maximum values of dose to assign to voxels with lowest value of cellularity and mode (most frequent) values of cellularity, respectively. These for us are 54 Gy (corresponding to dose assigned to CTV_low dose) and MAX_DOSE, specific to the patient.
        pdf_vals = tuple(zip(kde(X), X))
        N_dmax = max(pdf_vals)[1]
        dose = transfer_func(N0, Dmin, Dmax, N_dmax)
        
        #%%Plot the PDF and derived dose values as function of cellularity
        fig, ax = plt.subplots()
        plt.scatter(N0, dose, marker='x')
        ax2 = ax.twinx()
        ax.plot([N_dmax,N_dmax],[Dmin, max(dose)], 'k--')
        ax.plot([0,max(N0)],[Dmax, Dmax], 'k--')
        ax.annotate('Max Dose {0}'.format(Dmax), (max(N0)*0.7, Dmax*1.02))
        ax.annotate('mode of N density {0}'.format(int(N_dmax)), (N_dmax*1.1, Dmin*1.02))
        ax2.fill_between(X, PDF, fc="#AAAAFF", alpha=0.5)
        ax.set_ylabel('Dose [Gy]')
        ax2.set_ylabel('Probability density function')
        ax.set_xlabel('N of cells/voxel')
        plt.savefig(out_dir+'/PDF_DOSE_N_'+cell_type+'_plot.jpg')
        plt.close()
        
        #%% Compare new TCP with TCP that would be obtained by deliverying an homogeneous dose boost to the GTV to maximum dose estimated in the optimization
        TCP_heter =np.prod( np.exp( - N0 * np.exp(-alpha*dose-(beta*(dose**2))/n))) #TCP obtained with optimised heterogeneous dose distribution
        TCP_boost = np.prod( np.exp( - N0 * np.exp(-alpha*MAX_DOSE-(beta*(MAX_DOSE**2))/n))) # TCP with optimised homogeneous max dose boost to entire GTV
        print('TCP optimised heterogeneous dose plan: ',TCP_heter)
        print('TCP optimised homogeneous max dose boost plan: ',TCP_boost)
        
        plt.hist(dose, bins=20)
        plt.xlabel("Dose (Gy)")
        plt.ylabel("Counts")
        plt.savefig(out_dir+'/Dose_histo_'+cell_type+'.jpg')
        plt.close()
        
    #%% Saving the optimised dose values into a prescription dose image
        
        dose_nda = np.ndarray(shape=orig_shape) #Generate a 3D array of the same shape of the cell density map
        
        for d,i in zip(list(dose),listOfCoordinates): #Assign the calculated optimised dose to each voxel, based on the index of the corresponding cellularity voxel
            dose_nda[i]=d
    
        dose_img = sitk.GetImageFromArray(dose_nda) #Convert the numpy array into an image 
        dose_img.CopyInformation(cell_map) #Copy the headers of the cellularity image onto the headers of the optimised dose image
        
        sitk.WriteImage(dose_img, out_dir+'/Dose_optimised_CTV'+cell_type+'.nii') #Save the prescribed optimised dose image as a nifti file 
        
    #%%Set values within CTV with cellularity untrusted and bone pixels to receive Dmin = 54 Gy
        
        CTV_complete = sitk.ReadImage(data_supradir+current+'/MRI/baseline/orig/CTV_inDWI_ITK.nii')
        CTV_remaining = (CTV_complete>0) & (dose_img==0)
        dose_img_final = set_mask_value(dose_img, CTV_remaining, Dmin)
        sitk.WriteImage(dose_img_final, out_dir+'/Dose_optimised_CTV'+cell_type+'_final.nii')
        
    #%% Need to implement version where outliers vox in margins CTV receive 54 Gy and outtliers vox in CTV_HD receive 72 Gy
        
