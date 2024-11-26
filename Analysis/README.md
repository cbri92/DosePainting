## SBC Proton Dose Painting - Patients data analysis

### Author: Caterina Brighi

Analysis of the SBC patients radiotherapy planning data.
For a correct analysis pipeline, run scripts following the tasks order Task 1 - Task 5. 

#### Task 1: Convert RT dose plans from dicom to nifti: Conversion_RTplans_DICOM_to_nii.py

*    [x]  Convert dose painting and baseline RT dose dicom files into nifti format
         
#### Task 2: Plans conformity evaluation: Calculate QF and TCP.py

*    [x]  Calculate Quality Factor (QF) in GTV and CTV for both baseline and dose painting plans, and generate QF maps within the target volumes
*    [x]  Calculate Tumour Control Probability (TCP) in GTV and CTV for both baseline and dose painting plans, using dose plans and cellularity data from DWI MRI 

#### Task 3: Calculate dose statistics within the radiotherapy targets: CTV_GTV_doseStats.py

*    [x]  Calculate mean, std, median, minimum, maximum, D95% and D1% dose within the GTV and CTV for both baseline and dose painting plans

#### Task 4: Calculate DVH in radiotherapy targets and organs at risk: nii_RT_to_DVH.py

*    [x]  Calculate the DVH in radiotherapy targets and organs at risk for both baseline and dose painting plans

#### Task 5: Generate mean DVH plots for the entire patients' cohort in radiotherapy targets and organs at risk: Mean_DVH_plot.py

*    [x]  Generate mean DVH plots for both baseline and dose painting plans for the entire patients' cohort in radiotherapy targets and organs at risk. Also plot 25th and 75th percentile lines for percentage Volume. 
