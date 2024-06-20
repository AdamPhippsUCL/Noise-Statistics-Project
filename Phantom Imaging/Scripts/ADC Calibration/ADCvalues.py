
import numpy as np
import os
import sys
import SimpleITK as sitk


# Load ADC image
ADC = sitk.GetArrayFromImage(sitk.ReadImage(r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Outputs\20230601_QIBAPHANTOM\Images\dADC_corrected.mha"))

# Define ROI folder
ROIfolder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Outputs\20230601_QIBAPHANTOM\ROIs"

# Initialise Array for ADC results
ADCresults = np.zeros((8,3))

for ROInum in range(1,9):

    # Load ROI
    ROI = (sitk.GetArrayFromImage(sitk.ReadImage(f'{ROIfolder}/ROI{ROInum}.img')) ==1)

    ADCvals = ADC[ROI]

    ADCresults[ROInum-1, 0] = np.mean(ADCvals)
    ADCresults[ROInum-1,1] = np.percentile(ADCvals, 5)
    ADCresults[ROInum-1, 2] = np.percentile(ADCvals, 95)

print(ADCresults)