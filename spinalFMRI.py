def preInitialize():
    import nilearn

def inputs():
    fmri_img = 'E:\\Alex\\fMRI-spinal\\gul\\gul_WIP_FE_EPI_neck_20240612161710_501_t597000_moco.nii.gz'
    mask_img = 'E:\\Alex\\fMRI-spinal\\gul\\gul_WIP_3D_Spine_VIEW_T2W_20240612161710_401_seg_reg.nii.gz'
    labels_img = 'E:\\Alex\\fMRI-spinal\\gul\\gul_WIP_3D_Spine_VIEW_T2W_20240612161710_401_seg_labeled_discs.nii.gz'
    return fmri_img, mask_img, labels_img

def openImage(imagePath):
    from nilearn import image
    data = image.get_data(imagePath)
    return data 

def quantifyBOLD():
    #This script is dedicated to quantify signal in spinal fMRI images;
    import numpy as np

    myDataset = inputs()
    img = []
    for i in myDataset:
        img.append(openImage(i))
    
    masked_image = np.multiply(img[0], img[1])
    return img

img = quantifyBOLD()