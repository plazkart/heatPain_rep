function [BOLD_beta] = checkBOLDactivation(filepath, maskPath, metrics)
%CHECKBOLDACTIVATION Summary of this function goes here
% use as [~, metric] = hp_make('count_fmriMetrics', id, metrics)
% 'beta' - quantifying mean ratio of beta coef. 1 (stim) and
% constant term (last one)

%Check this for Visual stimulation data
% metrics = varargin{2};
% [filepath,name] = fileparts(filepath);

% callSmothing([filepath '\w' name '.nii']);
% V_mask = spm_vol([filepath '\sw' name '.nii']);

switch metrics
    case 'beta'
        fils = dir([mainStruct.meta.folder '\' nam '\derived\res\beta_*']);
        img_beta1= [fils(1).folder '\' fils(1).name];
        img_beta2 = [fils(end).folder '\' fils(end).name];
        V_img_1 = spm_vol(img_beta1);
        makeSameResolution(V_img_1.fname, V_mask.fname);
        img_mask = [filepath '\rsw' name '.nii'];

        BOLD_beta = countBeta(img_mask, img_beta1, img_beta2);
        mainStruct.(nam).proc.bold.mean_delta_beta1 = BOLD_beta;
        varargout{1} = BOLD_beta;
        hp_make('save', mainStruct);

    case 'insula_beta'
        fils = dir([filepath '\beta_*']);
%         img_mask = 'G:\_other\fMRS-heatPain\_meta\atlas_map\rsinsula_atlas.nii';
        img_mask = maskPath;  
        img_beta1= [fils(1).folder '\' fils(1).name];
%            makeSameResolution(img_beta1, img_mask)
        img_beta2 = [fils(end).folder '\' fils(end).name];
%         img_beta2 = [fils(end).folder '\' fils(18).name];
        BOLD_beta = countBeta(img_mask, img_beta1, img_beta2);
        varargout{1} = BOLD_beta;

    case 'cluster'
        fils = dir([mainStruct.meta.folder '\' nam '\derived\res\beta_*']);
        img_mask = [fils(1).folder '\clus1.nii'];
        img_beta1= [fils(1).folder '\' fils(1).name];
        img_beta2 = [fils(end).folder '\' fils(end).name];
        BOLD_beta = countBeta(img_mask, img_beta1, img_beta2);
        varargout{1} = BOLD_beta;


end
end


function BOLD_beta = countBeta(img_mask, img_beta1, img_beta2)
    %get filenames of the needed images
    V_img_1 = spm_vol(img_beta1);
    V_img_end = spm_vol(img_beta2);
    V_mask = spm_vol(img_mask);

    I_mask = V_mask.private.dat(:,:,:);
    I_img_1 = V_img_1.private.dat(:,:,:);
    I_img_end = V_img_end.private.dat(:,:,:);
    
    I_rel = I_img_1./I_img_end;
    MultiplyI = I_mask.*I_rel; idx = ~isnan(MultiplyI);
    MultiplyI = MultiplyI(idx); I_mask = I_mask(idx);
    BOLD_beta = sum(MultiplyI,'all')/sum(I_mask, 'all');

end

function makeSameResolution(img1, img2)
    nrun = 1;
    jobfile = {'C:\Users\Science\Documents\GitHub\matlab_work\My_scripts\utils\fmri_CoregReslice_job.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(2, nrun);
    for crun = 1:nrun
        inputs{1, crun} = {img1}; % Coregister: Reslice: Image Defining Space - cfg_files
        inputs{2, crun} = {img2}; % Coregister: Reslice: Images to Reslice - cfg_files
    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
end

