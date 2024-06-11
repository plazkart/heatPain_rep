function [mainStruct, varargout] = hpfmri_batch(action, varargin)
%HPFMRI_BATCH Summary of this function goes here
%   Sequence of analysis
% hpfmri_batch('new', mainStruct);
% hpfmri_batch('parseMRI', mainStruct, id);
% hpfmri_batch('parseFuncTable', mainStruct, id);

switch action
    case 'init'
        %use this just one time
        mainStruct = struct();
        mainStruct.meta.folder = 'G:\_other\fMRI-thermal\fmri_thermal';

        mainStruct.meta.subNumbers = 0;
        mainStruct.meta.date = datetime("today");
        save([mainStruct.meta.folder '\_meta\structMeta.mat'], 'mainStruct');

        
        %add new participant
    case 'new'
        if length(varargin)>0
            mainStruct = varargin{1};
            who_id = mainStruct.meta.subNumbers+1;
        end
        nam = sprintf('sub_%02i', who_id);
        mainStruct.(nam).id = mainStruct.meta.subNumbers+1;
        mainStruct.(nam).nam = [];
        mainStruct.(nam).date = datetime('today');
        mainStruct.(nam).data_check.sp = 0; %0 - no info, 1 - found, 2 - not found
        mainStruct.(nam).data_check.t1 = 0;
        mainStruct.(nam).data_check.fmri = 0;
        mainStruct.(nam).data_check.t2check = 0;
        mainStruct.(nam).folder = ['\' nam];

        mainStruct.(nam).proc.start_dynamic = 0;
        mainStruct.(nam).proc.tp_matrix = [];
        mainStruct.(nam).proc_check.tp_matrix = 0;
        mainStruct.(nam).proc.dummy_time = 0;
        
        mainStruct.(nam).proc_check.timepoints_spectra.sham = 0;
        mainStruct.(nam).proc_check.timepoints_spectra.act = 0;
        mainStruct.(nam).proc_check.tp_spectra_res.sham = 0;
        mainStruct.(nam).proc_check.tp_spectra_res.act = 0;
        
        mainStruct.meta.subNumbers = mainStruct.meta.subNumbers+1;

        mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\meta']);
        txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'w');
        fprintf(txt_protocol, 'Data processing steps log-file. \n Subject: %s \n Date: %s \n', nam, datetime("today"));
        fclose(txt_protocol);

        hpfmri_batch('save', mainStruct)

        % parsing data and place it into correct directory
    case 'save'
        %use as [mainStruct, varargout] = hpfmri_batch('save', mainStruct)
        if length(varargin)>0
            mainStruct = varargin{1};
        end
        mainStruct.meta.date = datetime("today");
        save([mainStruct.meta.folder '\_meta\structMeta.mat'], 'mainStruct');


    case 'load'
        default_HPfmri;

    case 'parseMRI'
        %use as [mainStruct, varargout] = hpfmri_batch('parseMRI', mainStruct, id)
        if length(varargin)>0
            mainStruct = varargin{1};
            if length(varargin)>1
                who_id = varargin{2};
            else
                who_id = mainStruct.meta.subNumbers;
            end
        end
        nam = sprintf('sub_%02i', who_id);
        fils_mri = dir([mainStruct.meta.folder mainStruct.(nam).folder '\*.nii']);
        ses = 1;
        if ~isempty(fils_mri)
            mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\func'])
            mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\anat'])
            
            for i=1:length(fils_mri)
                modalityI = regexp( fils_mri(i).name ,'T1W','match');
                if ~isempty(modalityI)
                    copyfile([mainStruct.meta.folder mainStruct.(nam).folder '\' fils_mri(i).name],...
                        [mainStruct.meta.folder mainStruct.(nam).folder '\anat\' nam '_anat.nii']);
                    mainStruct.(nam).data_check.t1 = 1;
                end
                modalityI = regexp( fils_mri(i).name ,'FE_EPI','match');
                if ~isempty(modalityI)
                    if mainStruct.(nam).data_check.fmri>0
                        mainStruct.(nam).data_check.fmri = mainStruct.(nam).data_check.fmri + 1;
                        ses = ses+1;
                    else
                        mainStruct.(nam).data_check.fmri = 1;
                    end
                    filenameFmri = sprintf('%s_func_%02i.nii', nam, ses);
                    copyfile([mainStruct.meta.folder mainStruct.(nam).folder '\' fils_mri(i).name],...
                        [mainStruct.meta.folder mainStruct.(nam).folder '\func\' filenameFmri ]);
                    mainStruct.(nam).proc.fmri{ses} = [mainStruct.meta.folder mainStruct.(nam).folder '\func\' filenameFmri];
                end
                
            end
            mainStruct.(nam).data_check.sesions = ses;
        end
        hpfmri_batch('save', mainStruct);

    case 'parseFuncTable'
        %use as mainStruct = hpfmri_batch('parseFuncTable', mainStruct, id);
        mainStruct = hpfmri_batch('load');
        if length(varargin)>0
%                         mainStruct = varargin{1};
            if length(varargin)>1
                who_id = varargin{2};
            else
                who_id = mainStruct.meta.subNumbers;
            end
        end
        nam = sprintf('sub_%02i', who_id);
        fils_xl = dir([mainStruct.meta.folder mainStruct.(nam).folder '\*.xlsx']);
        if length(fils_xl)>1
            mainStruct.(nam).funcTable.pretest = fils_xl(1).name;
            mainStruct.(nam).data_check.funcTable = 1;
            mainStruct.(nam).funcTable.fmri = fils_xl(2).name;
            mainStruct.(nam).data_check.funcTable =mainStruct.(nam).data_check.funcTable+10;
            ses_nam = sprintf('_ses_%02i', 1);
            copyfile([fils_xl(2).folder '\' fils_xl(2).name], [mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold' ses_nam '.xlsx']);
            if length(fils_xl)>2
                mainStruct.(nam).funcTable.fmrs = fils_xl(3).name;
                mainStruct.(nam).data_check.funcTable =mainStruct.(nam).data_check.funcTable+100;
                ses_nam = sprintf('_ses_%02i', 2);
                copyfile([fils_xl(3).folder '\' fils_xl(3).name], [mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold' ses_nam '.xlsx']);
            else
                mainStruct.(nam).funcTable.fmrs = 'n';
            end
        else
            mainStruct.(nam).data_check.funcTable =0;
        end
        fils_csv = dir([mainStruct.meta.folder mainStruct.(nam).folder '\*.csv']);
        if length(fils_csv)>1
            mainStruct.(nam).funcTable.est = fils_csv(1).name;
            copyfile([fils_csv(1).folder '\' fils_csv(1).name], [mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_est.csv']);
        end
        hpfmri_batch('save', mainStruct);


        case 'fmri_proc'
        % use it as hpfmri_batch('fmri_proc', id, mainStruct, procSteps); mainStruct - if
        % needed
        if length(varargin)>0
            id = varargin{1};
            if length(varargin)>1
                mainStruct = varargin{2};
            else
                mainStruct = hp_make('load');
            end
        else
            error('There no subject name to process');
        end
        nam = sprintf('sub_%02i', id);

        % rewrite as needed 

        if mainStruct.(nam).data_check.fmri > 0
            for i=1:length(mainStruct.(nam).proc.fmri)
                fmri_img = spm_vol(mainStruct.(nam).proc.fmri{i});
                NSA(i) = length(fmri_img);
                for ii=1:NSA
                    func_data{1, i}{ii, 1} = [mainStruct.(nam).proc.fmri{i} ',' num2str(ii)];
                end
            end
        end
        if mainStruct.(nam).data_check.t1 == 1
            anat_data(1, 1) = {[mainStruct.meta.folder mainStruct.(nam).folder '\anat\' nam '_anat.nii']};
        end
        
        spatial_steps = 1;
        if length(varargin)>2
            %varargin{3}: if 0, not to do spatial steps, if 1 -yes
            if varargin{3}<1
                spatial_steps = 0;
            end
        end
        
        if spatial_steps
            nrun = 1;
%             jobfile = {[mainStruct.meta.folder '\_meta\spatial_empty_job.m']};
%             jobs = repmat(jobfile, 1, nrun);
            
            inputs = cell(3, nrun);
            for crun = 1:nrun
                inputs{1, crun} = func_data; % Realign: Estimate & Reslice: Session - cfg_files
                inputs{2, crun} = anat_data; % Coregister: Estimate: Reference Image - cfg_files
                inputs{3, crun} = anat_data; % Segment: Volumes - cfg_files
            end
            
%             spm('defaults', 'FMRI');
%             spm_jobman('run', jobs, inputs{:});
             callfMRIProcessing(inputs, 'interleaved', 1);

            %There are new files in the fmri-directory, so we will copy
            %them into another dir (called derived). Initial files are
            %saved
            fils_func = dir([mainStruct.meta.folder mainStruct.(nam).folder '\func\*' nam '*']);
            mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\derived']);
            for i=1:length(fils_func)
                dontdelete_flag = 0;
                if ~contains(fils_func(i).name, 'xlsx')
                    for ii =1:length(mainStruct.(nam).proc.fmri)
                        [~, fl_name] = fileparts(mainStruct.(nam).proc.fmri{ii});
                        if strcmp(fils_func(i).name, [fl_name '.nii']) || dontdelete_flag==1
                            dontdelete_flag = 1;
                        end
                    end
                    if dontdelete_flag == 0
                        copyfile([fils_func(i).folder '\' fils_func(i).name], ...
                            [mainStruct.meta.folder mainStruct.(nam).folder '\derived\' fils_func(i).name]);
                        delete([fils_func(i).folder '\' fils_func(i).name]);
                    end
                end
            end
            mainStruct.(nam).data_check.fmri_spat = 1;
        end

        mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\derived\res']);
        nrun = 1; % enter the number of runs here
        jobfile = {[mainStruct.meta.folder '\_meta\stats_empty_job.m']};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(5, nrun);
        for crun = 1:nrun
            inputs{1, crun} = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\res']}; % fMRI model specification: Directory - cfg_files
            for ii =1:length(mainStruct.(nam).proc.fmri)
                func_data = cell(1);
                ses_nam = sprintf('_ses_%02i', ii);
                [~, fl_name] = fileparts(mainStruct.(nam).proc.fmri{ii});
                for i=1:NSA(ii)
                    func_data{1, i} = [mainStruct.meta.folder mainStruct.(nam).folder '\derived\swra' fl_name '.nii,' num2str(i)];
                end
                inputs{2, ii} = func_data; % fMRI model specification: Scans - cfg_files
                if floor(mod(mainStruct.(nam).data_check.funcTable, 100)/10)>0
                    regressorList = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold' ses_nam '.xlsx']);
                    startTime = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold' ses_nam '.xlsx']);
                     regressorList(:, 1) = regressorList(:, 1) - startTime(1) - (12 - mainStruct.(nam).proc.dummy_time);
%                     regressorList(:, 1) = regressorList(:, 1) - startTime(1);
                end
                inputs{3, ii} = regressorList(:, 1); % fMRI model specification: Onsets - cfg_entry
                inputs{4, ii} = regressorList(:, 2); % fMRI model specification: Durations - cfg_entry
                inputs{5, ii} = [mainStruct.meta.folder mainStruct.(nam).folder '\derived\rp_a' fl_name '.txt']; % fMRI model specification: Multiple regressors - cfg_files
            end
        end
        callGLM(inputs);
        

end
end

function callfMRIProcessing(inputs, slice_case, segment)
    mainStruct = hpfmri_batch('load');
    
    k=1;
    sess = length(inputs{1, 1});
    an_image = spm_vol(inputs{1, 1}{1}{1});
    img_info = niftiinfo(an_image.fname);

    matlabbatch{k}.spm.temporal.st.scans = {inputs{1, 1}{:}};
    matlabbatch{k}.spm.temporal.st.nslices = img_info.ImageSize(3);
    matlabbatch{k}.spm.temporal.st.tr = img_info.PixelDimensions(4);
    matlabbatch{k}.spm.temporal.st.ta = matlabbatch{k}.spm.temporal.st.tr - matlabbatch{k}.spm.temporal.st.tr/matlabbatch{k}.spm.temporal.st.nslices;
    switch slice_case
        case 'interleaved'
            matlabbatch{1}.spm.temporal.st.so = [1	7	13	19	25	31	2	8	14	20	26	32	3	9	15	21	27	33	4	10	16	22	28	34	5	11	17	23	29	35	6	12	18	24	30];
        case 'ascending'
            matlabbatch{1}.spm.temporal.st.so = [1:matlabbatch{k}.spm.temporal.st.nslices];
    end
    matlabbatch{k}.spm.temporal.st.refslice = 1;
    matlabbatch{k}.spm.temporal.st.prefix = 'a';
    k=k+1;
    matlabbatch{k}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    if sess>1
        matlabbatch{k}.spm.spatial.realign.estwrite.data{2}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 2)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
    end
    matlabbatch{k}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{k}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{k}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{k}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{k}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{k}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{k}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{k}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{k}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{k}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{k}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{k}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    k=k+1;
    matlabbatch{k}.spm.spatial.coreg.estimate.ref = inputs{2, 1};
    matlabbatch{k}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
%     matlabbatch{3}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
% matlabbatch{3}.spm.spatial.coreg.estimate.other(2) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
    matlabbatch{k}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
    if sess>1
        matlabbatch{k}.spm.spatial.coreg.estimate.other(2) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
    end
    
    matlabbatch{k}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{k}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{k}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{k}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    if segment
        k=k+1;
        matlabbatch{k}.spm.spatial.preproc.channel.vols = inputs{3, 1};
        matlabbatch{k}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{k}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{k}.spm.spatial.preproc.channel.write = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(1).tpm = {[mainStruct.meta.SPMfolder '\spm12\tpm\TPM.nii,1']};
        matlabbatch{k}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{k}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{k}.spm.spatial.preproc.tissue(1).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(2).tpm = {[mainStruct.meta.SPMfolder '\spm12\tpm\TPM.nii,2']};
        matlabbatch{k}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{k}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{k}.spm.spatial.preproc.tissue(2).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(3).tpm = {[mainStruct.meta.SPMfolder '\spm12\tpm\TPM.nii,3']};
        matlabbatch{k}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{k}.spm.spatial.preproc.tissue(3).warped = [1 1];
        matlabbatch{k}.spm.spatial.preproc.tissue(4).tpm = {[mainStruct.meta.SPMfolder '\spm12\tpm\TPM.nii,4']};
        matlabbatch{k}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{k}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{k}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{k}.spm.spatial.preproc.tissue(5).tpm = {[mainStruct.meta.SPMfolder '\spm12\tpm\TPM.nii,5']};
        matlabbatch{k}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{k}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{k}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{k}.spm.spatial.preproc.tissue(6).tpm = {[mainStruct.meta.SPMfolder '\spm12\tpm\TPM.nii,6']};
        matlabbatch{k}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{k}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{k}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{k}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{k}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{k}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{k}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{k}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{k}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{k}.spm.spatial.preproc.warp.write = [0 1];
        matlabbatch{k}.spm.spatial.preproc.warp.vox = NaN;
        matlabbatch{k}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
            NaN NaN NaN];
    end
    k=k+1;
    if segment
        matlabbatch{k}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    else
        temp_dir = split(inputs{2}, '\');
        temp_dir{end} = ['y_' temp_dir{end}];
        new_name = fullfile(temp_dir{:});
        matlabbatch{k}.spm.spatial.normalise.write.subj.def(1) = {new_name};
    end
    matlabbatch{k}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{k}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{k}.spm.spatial.normalise.write.woptions.vox = [1.43 1.43 3];
    matlabbatch{k}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{k}.spm.spatial.normalise.write.woptions.prefix = 'w';
    k=k+1;
    if segment
        matlabbatch{k}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        else
        matlabbatch{5}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    end
    matlabbatch{k}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{k}.spm.spatial.smooth.dtype = 0;
    matlabbatch{k}.spm.spatial.smooth.im = 0;
    matlabbatch{k}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('run',matlabbatch);

end

function callGLM(inputs)
    matlabbatch{1}.spm.stats.fmri_spec.dir = inputs{1, 1};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
for ii=1:size(inputs, 2)
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).scans = inputs{2, ii}';
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).cond.name = 'thermal';
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).cond.onset = inputs{3, ii};
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).cond.duration = inputs{4, ii};
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).cond.tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).cond.orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).multi_reg = {inputs{5, ii}};
    matlabbatch{1}.spm.stats.fmri_spec.sess(ii).hpf = 128;
end
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'thermal';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'both';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.ps = true;

spm_jobman('run',matlabbatch);
end


function startTime = getTTLtime(xlsfile)
% stimulus times ONLY for MRS modality
    k=1;
    eventsTable = readtable(xlsfile);
    if iscell(eventsTable.Events)
        TTL = strfind(eventsTable.Events, 'TTL');
        for i=1:length(TTL)
            if ~isempty(TTL{i})
                startTime(k) = eventsTable.Timestamp_msec_(i)/1000;
                k=k+1;
%                 break
            end
        end
    else
        tstamp_1 = eventsTable.Tec_C_-35;
        thrs_1 = find(tstamp_1>0.1); thrs_1 = thrs_1(1);
        startTime = eventsTable.Timestamp_msec_(thrs_1)/1000;
    end
    

end