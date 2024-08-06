def preInitialize():
    import nilearn

def inputs():
    fmri_img = 'E:\\Alex\\fMRI-spinal\\gul\\gul_WIP_FE_EPI_neck_20240612161710_501_t597000_moco.nii.gz'
    mask_img = 'E:\\Alex\\fMRI-spinal\\gul\\gul_WIP_3D_Spine_VIEW_T2W_20240612161710_401_seg_reg.nii.gz'
    labels_img = 'E:\\Alex\\fMRI-spinal\\gul\\gul_WIP_3D_Spine_VIEW_T2W_20240612161710_401_seg_labeled.nii.gz'
    return fmri_img, mask_img, labels_img

def openImage(imagePath):
    from nilearn import image
    data = image.load_img(imagePath)
    return data 

def applyMask(img, mask_img):

    from nilearn.masking import apply_mask
    from nilearn import image

    masked_data = apply_mask(img, mask_img)

    return masked_data

    # masked_data shape is (timepoints, voxels). We can plot the first 150
    # timepoints from two voxels
        

    # sphinx_gallery_dummy_images=3

def quantifyBOLD():
    #This script is dedicated to quantify signal in spinal fMRI images;
    import numpy as np
    from nilearn import image
    from nilearn import masking

    myDataset = inputs()
    img = []
    for i in myDataset:
        img.append(openImage(i))
    #resample image and get one vertebrae
    img[1] = image.binarize_img(img[1],  threshold=0.99)
    masked_image = []
    for i in range(5, 9):
        img_label = image.math_img("img == "+str(i), img=img[2])
        img_label = image.resample_to_img(img_label, img[1], interpolation= "nearest")

        #intersect two masks
        mask_img = masking.intersect_masks([img_label, img[1]], threshold=0.5, connected=True)
        masked_image.append(applyMask(img[0], mask_img))

    return masked_image

def plotTheShit(timecurve):
        import numpy as np
        v_x = np.zeros([200, len(timecurve)])
        k = 0
        for i in timecurve:
             a = i.mean(axis = 1)-i.mean()
             a = a/i.shape[1]
             v_x[:, k] = a
             k = k + 1
             
    # And now plot a few of these
        import matplotlib.pyplot as plt

        plt.figure(figsize=(7, 5))
        for idx in range(4):
            plt.plot(v_x[:, idx])
        plt.xlabel("Time [TRs]", fontsize=16)
        plt.ylabel("Intensity", fontsize=16)
        plt.xlim(0, 150)
        plt.subplots_adjust(bottom=0.12, top=0.95, right=0.95, left=0.12)

        plt.show()

timecurve = quantifyBOLD()
plotTheShit(timecurve)