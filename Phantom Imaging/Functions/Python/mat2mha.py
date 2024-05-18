from scipy.io import loadmat
import SimpleITK as sitk
import numpy as np
import os


# MAT filename
Matfname = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\New VERDICT outputs\Original VERDICT\Original ex905003000\AMICO\shV_20240307\fIC.mat"

# Get folder path and img name
folder, imgname = os.path.split(Matfname) 
imgname = imgname[0:-4] # Removing .mat

# Get image array
img = (loadmat(Matfname))[imgname]

# remove infinities
img[img == np.inf] = 0
# Remove nan
img[np.isnan(img)] = 0

# Change image orientation
img = np.moveaxis( img , -1, 0)  



# Save as mha
sitk.WriteImage(sitk.GetImageFromArray(img), f'{folder}/{imgname}.mha')
