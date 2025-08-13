def preInitialize():
    import nilearn

class imagesFMRI:
    epi = ''
    spinalCord = ''
    labelsCord = ''
    label1 = ''
    label3 = ''
    label5 = ''
    label35 = ''
    label4 = ''
    label8 = ''
    label12 = ''
    label30 = ''

def inputs():
    #parse all needed data from the spinal fmri directory
    import os
    from nilearn import image
    
    mainDir = 'E:\\Alex\\fMRI-spinal\\Ivantsov\\spinal-fmri'
    filsList = os.listdir(mainDir)
    os.chdir(mainDir)
    dats = imagesFMRI()
    for i in filsList:
        if i.find('epi_moco.nii.gz')>-1:
            dats.epi = image.load_img(i)
        if i.find('t2_seg_reg.nii.gz')>-1:
            dats.spinalCord = image.load_img(i)
    filsList = os.listdir(mainDir + '\\label\\template')
    for i in filsList:
        if i.find('spinal_levels.nii.gz')>-1:
                dats.labelsCord = image.load_img('label\\template\\' + i)
    filsList = os.listdir(mainDir + '\\label\\atlas')
    atlasLabels = [1, 3, 5, 35, 4, 8, 12, 30]
    for i in filsList:
        for j in atlasLabels:
            if i.find('atlas_{:02}.nii.gz'.format(j))>-1:
                exec('dats.label' + str(j) + ' = image.load_img(\'label\\\\atlas\\' + i + '\')')
    
    return dats


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
    import os
    from nilearn import image
    from nilearn import masking

    myDataset = inputs()
    unifiedAtlases = list()
    #make mask of activated area
    unifiedAtlases.append(image.math_img("img1 + img2+ img3 + img4",
                      img1=myDataset.label1, img2=myDataset.label3,
                      img3=myDataset.label5, img4=myDataset.label35))
    unifiedAtlases[0] = image.binarize_img(unifiedAtlases[0],  threshold=0.01)

    unifiedAtlases.append( image.math_img("img1 + img2+ img3 + img4",
                      img1=myDataset.label4, img2=myDataset.label8,
                      img3=myDataset.label12, img4=myDataset.label30))
    unifiedAtlases[1] = image.binarize_img(unifiedAtlases[1],  threshold=0.01)

    masked_image = []
    for j in range(2):
        for i in range(6, 9):
            img_label = image.math_img("img == "+str(i), img=myDataset.labelsCord)
            img_mask = masking.intersect_masks([img_label, unifiedAtlases[j]], threshold=0.5, connected=True)
            img_mask = image.resample_to_img(img_mask, myDataset.epi, interpolation= "nearest")

            #intersect two masks
    #       mask_img = masking.intersect_masks([img_mask, myDataset.epi], threshold=0.5, connected=True)
            masked_image.append(applyMask(myDataset.epi, img_mask))

    return masked_image

class resultsExp:
    SNR = []

def qunatificationSteps(timecurve):
    import numpy as np
    out = resultsExp
    # count SNR in all regions
    for i in timecurve:
        SNR = np.mean(i)/np.std(i)
        out.SNR.append(SNR)
    #count 
    
    return out

def plotTheShit(timecurve):
        import numpy as np
        v_x = np.zeros([250, len(timecurve)])
        k = 0
        for i in timecurve:
             a = i.mean(axis = 1)-i.mean()
             a = a/i.shape[1]
             v_x[:, k] = a
             k = k + 1
             
    # And now plot a few of these
        import matplotlib.pyplot as plt

        plt.figure(figsize=(7, 5))
        for idx in range(len(timecurve)):
            plt.plot(v_x[:, idx])
        plt.xlabel("Time [TRs]", fontsize=16)
        plt.ylabel("Intensity", fontsize=16)
        plt.xlim(0, 250)
        plt.subplots_adjust(bottom=0.12, top=0.95, right=0.95, left=0.12)

        plt.show()

def processSCANPHYSLOGdata(scanphyslogName, outDirName, fmri_file):
    from scanphyslog2bids.core import PhilipsPhysioLog
    import nibabel as nib
    import numpy as np

    log_file = scanphyslogName
    out_dir = outDirName  # where the BIDSified data should be saved
    deriv_dir = outDirName + '/physio'  # where some QC plots should be saved

    # fmri_file is used to extract metadata, such as TR and number of volumes
    fmri_file = fmri_file 
    fmri_img = nib.load(fmri_file)
    n_dyns = fmri_img.shape[-1]
    tr = np.round(fmri_img.header['pixdim'][4], 3)

    # Create PhilipsPhysioLog object with info
    phlog = PhilipsPhysioLog(f=log_file, tr=tr, n_dyns=n_dyns, sf=496, manually_stopped=False)

    # Load in data, do some preprocessing
    phlog.load()

    # Try to align physio data with scan data, using a particular method
    # (either "vol_markers", "gradient_log", or "interpolation")
    phlog.align(trigger_method='interpolate', which_grad='y')  # load and find vol triggers

    # Write out BIDS files
    phlog.to_bids(out_dir)  # writes out .tsv.gz and .json files

    # Optional: plot some QC graphs for alignment and actual traces
    phlog.plot_alignment(out_dir=deriv_dir)  # plots alignment with gradient
    phlog.plot_traces(out_dir=deriv_dir)  # plots cardiac/resp traces

timecurve = quantifyBOLD()
plotTheShit(timecurve)
#res = qunatificationSteps(timecurve)
#print(res.SNR)

#processSCANPHYSLOGdata('E:\\Alex\\fMRI-spinal\\Ivantsov\\spinal-fmri\\physIO\\ScanPsaLog20240906180252.log',
#                       'E:\\Alex\\fMRI-spinal\\Ivantsov\\spinal-fmri\\physIO' , 'E:\\Alex\\fMRI-spinal\\Ivantsov\\spinal-fmri\\Ivantsov_WIP_FE_EPI_neck_20240906170706_901.nii.gz')