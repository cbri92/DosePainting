## SBC Proton Dose Painting - Patients data pre-processing

### Author: Caterina Brighi

Pre-processing of the SBC patients data.
For a correct preprocessing pipeline, run scripts following the tasks order Task 1 - Task 6. 

#### Task 1: Data conversion: dicom2nifti.py

*    [x]  Convert MRI, RTDOSE, RTSTRUCT and CT DICOM data to NIFTI format

#### Task 2: Coregister RT data to DWI MRI space: Transform_RBEinDWI.py, Transform_GTVinDWIspace.py and Transform_CTVinDWIspace.py

*    [x]  Register CT images, RBE dose and GTV and CTV to T2 MRI image space, and then to DWI image space, by means applying given transformation matrices

#### Task 3: Patient's specific dose painting prescriptions derivation: Dpresc_fromCell_SBCproton.py

*    [x]  Generate a prescription dose map in DWI space, by linearly scaling the dose between 74-81 Gy with per-voxel cellularity in the GTV

#### Task 4: Patient's specific inverse dose painting prescriptions derivation: Inverse_DosePrescriptions_SBCproton.py

*    [x]  Register the prescription dose map from DWI to CT space by applying the inverse transforms from Task 2
*    [x]  Generate inverse dose prescription maps in the GTV as Dinv = Dmax - Dpresc, where Dmax = maximum dose estimated in GTV for each patient
*    [x]  Derive CTV high dose (HD) margins and assign 74 Gy to all voxels in HD margins
*    [x]  Obtain a Dose prescription and an inverse dose prescription map in nifitn format
*    [x]  Convert the inverse dose prescription nifti file into a RTdose dicom file



