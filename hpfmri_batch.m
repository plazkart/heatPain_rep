function [mainStruct, varargout] = hpfmri_batch(action, varargin)
%HPFMRI_BATCH Summary of this function goes here
%   Sequence of analysis
% hpfmri_batch('new', mainStruct);
% hpfmri_batch('parseMRI', mainStruct, id);
% hpfmri_batch('parseFuncTable', mainStruct, id);

%new pipeline
% hpfmri_batch('newBIDS');
% hpfmri_batch('parseBIDS', id)
% hpfmri_batch('parseXLSX', id);
% hpfmri_batch('fmri_procBIDS', id, procSteps);

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
        %use as [mainStruct, varargout] = hpfmri_batch('parseMRI', id)
        mainStruct = hpfmri_batch('load');
        who_id = varargin{1};
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
        %use as mainStruct = hpfmri_batch('parseFuncTable', id);
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


       

        %% Hereby all scripts are mado for processing in BIDSy-style
        case 'initAsBIDS'
        %use this just one time
        mainStruct = struct();
        mainStruct.meta.folder = 'E:\Alex\fmri_thermal';
        mainStruct.meta.bids_folder = [mainStruct.meta.folder '\_bids'];

        mainStruct.meta.subNumbers = 0;
        mainStruct.meta.date = datetime("today");
        

        mkdir([mainStruct.meta.folder '\_bids']);
        mkdir([mainStruct.meta.folder '\_bids\sourcedir']);
        mkdir([mainStruct.meta.folder '\_bids\derivatives']);
        mkdir([mainStruct.meta.folder '\_bids\derivatives\spm']);
        
        fil = fopen([mainStruct.meta.folder '\_bids\.bidsignore'], 'w');
        fclose(fil);
        fil = fopen([mainStruct.meta.folder '\_bids\README'], 'w');
        fclose(fil);

        mainStruct.meta.JSONstruct.Name = 'fMRI of thermal pain';
        mainStruct.meta.JSONstruct.BIDSVersion = '';
        mainStruct.meta.JSONstruct.License = '';
        mainStruct.meta.JSONstruct.Authors = {'A. Yakovlev', 'M. Ublinskyi', 'O. Bozhko', 'T. Akhadov'};
        mainStruct.meta.JSONstruct.HowToAcknowledge = ''; 
        mainStruct.meta.JSONstruct.Funding = '';
        mainStruct.meta.JSONstruct.DatasetDOI = '';
        txt = jsonencode(mainStruct.meta.JSONstruct);
        fil = fopen([mainStruct.meta.folder '\_bids\dataset_description.json'], 'w');
        fprintf(fil, '%s', txt);
        fclose(fil);
        save([mainStruct.meta.folder '\_meta\structMetaBIDS.mat'], 'mainStruct');

        case 'saveBIDS'
            %use as [mainStruct, varargout] = hpfmri_batch('saveBIDS', mainStruct)
            if length(varargin)>0
                mainStruct = varargin{1};
            end
            mainStruct.meta.date = datetime("today");
            save([mainStruct.meta.folder '\_meta\structMetaBIDS.mat'], 'mainStruct');


        case 'loadBIDS'
            % use as mainStruct = hpfmri_batch('loadBIDS');
            default_HPfmri;
            load([mainStruct.meta.folder '\_meta\structMetaBIDS.mat']);
            varargout{1} = mainStruct;

        case 'newBIDS'
            % use as mainStruct = hpfmri_batch('newBIDS');
            mainStruct = hpfmri_batch('loadBIDS');
            who_id = mainStruct.meta.subNumbers+1;
        nam = sprintf('sub_%02i', who_id);
        BIDSnam = sprintf('sub-%02i', who_id);

        mainStruct.(nam).id = mainStruct.meta.subNumbers+1;
        mainStruct.(nam).nam = [];
        mainStruct.(nam).date = datetime('today');
        mainStruct.(nam).data_check.sp = 0; %0 - no info, 1 - found, 2 - not found
        mainStruct.(nam).data_check.t1 = 0;
        mainStruct.(nam).data_check.fmri = 0;
        mainStruct.(nam).data_check.t2check = 0;
        mainStruct.(nam).folder = ['\' BIDSnam];

        mainStruct.(nam).proc.start_dynamic = 0;
        mainStruct.(nam).proc.tp_matrix = [];
        mainStruct.(nam).proc_check.tp_matrix = 0;
        mainStruct.(nam).proc.dummy_time = 12;
        
        mainStruct.meta.subNumbers = mainStruct.meta.subNumbers+1;
        
        mkdir([mainStruct.meta.folder '\_bids\' mainStruct.(nam).folder]);
        mkdir([mainStruct.meta.folder '\_bids\sourcedir\' mainStruct.(nam).folder]);
%         mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\meta']);
%         txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'w');
%         fprintf(txt_protocol, 'Data processing steps log-file. \n Subject: %s \n Date: %s \n', nam, datetime("today"));
%         fclose(txt_protocol);

        hpfmri_batch('saveBIDS', mainStruct);

    case 'parseBIDS'
        %use as [mainStruct, varargout] = hpfmri_batch('parseBIDS', id)
        mainStruct = hpfmri_batch('loadBIDS');
        id = varargin{1};
        
        nam = sprintf('sub_%02i', id);
        BIDSnam = sprintf('sub-%02i', id);
        sourceDIR = [mainStruct.meta.folder '\_bids\sourcedir\'  BIDSnam];

        %looking for mri-files
        fils = dir([sourceDIR '\*.nii']);
        if ~isempty(fils)
            mkdir([mainStruct.meta.bids_folder mainStruct.(nam).folder '\anat'])
            mkdir([mainStruct.meta.bids_folder mainStruct.(nam).folder '\func'])
            %             mkdir([mainStruct.meta.folder mainStruct.meta.bids_folder mainStruct.(nam).folder '\derived'])
            %first get info about converted files then copy as needed
            mriFiles = struct('img', [], 'type', []);img_i = 1;
            for i=1:length(fils)
                modalityI = regexp( fils(i).name ,'T1W','match');
                if ~isempty(modalityI)
                    mriFiles(img_i).img = [fils(i).folder '\' fils(i).name];
                    mriFiles(img_i).type = 'T1w'; img_i = img_i + 1;
                end
                modalityI = regexp( fils(i).name ,'EPI','match');
                if ~isempty(modalityI)
                    mriFiles(img_i).img = [fils(i).folder '\' fils(i).name];
                    mriFiles(img_i).type = 'EPI'; img_i = img_i + 1; 
                end
            end
            runNum = 1;
            for i =1:length(mriFiles)
                if contains(mriFiles(i).type, 'T1w')
                    copyfile(mriFiles(i).img, [mainStruct.meta.bids_folder mainStruct.(nam).folder '\anat\' BIDSnam '_T1w.nii'])
                    copyfile([mriFiles(i).img(1:end-3) 'json' ], [mainStruct.meta.bids_folder mainStruct.(nam).folder '\anat\' BIDSnam '_T1w.json'])
                    mainStruct.(nam).data_check.t1 = 1;
                end
                if contains(mriFiles(i).type, 'EPI')
                    runNam = sprintf('run-%02i', runNum); 
                    copyfile(mriFiles(i).img, [mainStruct.meta.bids_folder mainStruct.(nam).folder '\func\' BIDSnam '_task-heatpain_' runNam '_bold.nii'])
                    copyfile([mriFiles(i).img(1:end-3) 'json' ], [mainStruct.meta.bids_folder mainStruct.(nam).folder '\func\' BIDSnam '_task-heatpain_' runNam '_bold.json'])
                    mainStruct.(nam).data_check.EPI = runNum; runNum = runNum +1;
                end
            end
        end
        hpfmri_batch('saveBIDS', mainStruct);

    case 'parseXLSX'
        %use as [mainStruct, varargout] = hpfmri_batch('parseXLSX', id)
         %looking for xlsx-files
         mainStruct = hpfmri_batch('loadBIDS');
        
        id = varargin{1};
        nam = sprintf('sub_%02i', id);
        BIDSnam = sprintf('sub-%02i', id);
        sourceDIR = [mainStruct.meta.folder '\_bids\sourcedir\'  BIDSnam];

        fils = dir([sourceDIR '\*.xlsx']);
        for i=1:length(fils)
            releasedTime(i, 1) = fils(i).datenum;
        end
        [~, Idxs] = sort(releasedTime);
        k = 1;
        for i=1:length(Idxs)
            if fils(Idxs(i)).bytes > 200000
                reggr = heatPain_makeRegressor([fils(Idxs(i)).folder '\' fils(Idxs(i)).name]);
                startTime = getTTLtime([fils(Idxs(i)).folder '\' fils(Idxs(i)).name]);
                reggr(:, 1) = reggr(:, 1) - startTime(1);
                runNam = sprintf('run-%02i', k); 
                reggr = array2table(reggr,"VariableNames",{'onset', 'duration'});
                writetable(reggr, [mainStruct.meta.bids_folder mainStruct.(nam).folder '\func\' BIDSnam '_task-heatpain_' runNam '_events.tsv'],...
                    'Delimiter','tab', 'FileType','text');
                k = k+1;
            end
        end

         case 'fmri_procBIDS'
        % use it as hpfmri_batch('fmri_procBIDS', id, procSteps); mainStruct - if
        % needed
        mainStruct = hpfmri_batch('loadBIDS');
        if length(varargin)>0
            id = varargin{1};
        else
            error('There no subject name to process');
        end
        nam = sprintf('sub_%02i', id);
        BIDSnam = sprintf('sub-%02i', id);
        sesNum = mainStruct.(nam).data_check.EPI;

        %create folder to store processed data
        mkdir([mainStruct.meta.bids_folder 'derivatives\spm\' BIDSnam]);
        derivativesDir = [mainStruct.meta.bids_folder 'derivatives\spm\' BIDSnam];

        % rewrite as needed 
        spatial_steps = 1;
        if length(varargin)>1
            %varargin{3}: if 0, not to do spatial steps, if 1 -yes
            if varargin{2}<1
                spatial_steps = 0;
            end
        end

        %set inputs
        for i=1:sesNum
            runNam = sprintf('run-%02i', i);
            fmri_name = [mainStruct.meta.bids_folder BIDSnam '\func\' BIDSnam '_task-heatpain_' runNam '_bold.nii'];
            fmri_img = spm_vol(fmri_name);
            NSA(i) = length(fmri_img);
            for ii=1:NSA
                func_data{1, i}{ii, 1} = [fmri_name ',' num2str(ii)];
            end
        end
        anat_data(1, 1) = {[mainStruct.meta.bids_folder BIDSnam  '\anat\' BIDSnam '_T1w.nii']};
        %save initial filenames
        fl_name = dir([mainStruct.meta.bids_folder BIDSnam '\*\' BIDSnam '*']);

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
%             callfMRIProcessing(inputs, 'ascending', 1); 

            %There are new files in the fmri-directory, so we will copy
            %them into another dir (called derived). Initial files are
            %saved
            fils_func = dir([mainStruct.meta.bids_folder mainStruct.(nam).folder '\*\*' BIDSnam '*']);

            for i=1:length(fils_func)
                dontdelete_flag = 0;
                for ii = 1:length(fl_name)
                    [~, tempName, ext] = fileparts(fl_name(ii).name);
                    if strcmp(fils_func(i).name, [tempName, ext])
                        dontdelete_flag = 1;
                    end
                end
                if dontdelete_flag == 0
                    copyfile([fils_func(i).folder '\' fils_func(i).name], ...
                        [derivativesDir '\' fils_func(i).name]);
                    delete([fils_func(i).folder '\' fils_func(i).name]);
                end
            end
            mainStruct.(nam).data_check.fmri_spat = 1;
        end

        mkdir([derivativesDir '\res']);
        nrun = 1; % enter the number of runs here
        jobfile = {[mainStruct.meta.folder '\_meta\stats_empty_job.m']};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(5, nrun);
        for crun = 1:nrun
            inputs{1, crun} = {[derivativesDir '\res']}; % fMRI model specification: Directory - cfg_files
            for ii =1:sesNum
                func_data = cell(1);
                runNam = sprintf('run-%02i', ii);
                fmri_name = [mainStruct.meta.bids_folder mainStruct.(nam).folder '\func\' BIDSnam '_task-heatpain_' runNam '_bold.nii'];

                [~, temp_name] = fileparts(fmri_name);
                for i=1:NSA(ii)
                    func_data{1, i} = [derivativesDir '\swra' temp_name '.nii,' num2str(i)];
                end
                inputs{2, ii} = func_data; % fMRI model specification: Scans - cfg_files
                regrr = readtable([mainStruct.meta.bids_folder mainStruct.(nam).folder '\func\' BIDSnam '_task-heatpain_' runNam '_events.tsv'], 'Delimiter','tab', 'FileType','text');
                regrr = table2array(regrr);
                inputs{3, ii} = regrr(:, 1); % fMRI model specification: Onsets - cfg_entry
                inputs{4, ii} = regrr(:, 2); % fMRI model specification: Durations - cfg_entry
                inputs{5, ii} = [derivativesDir '\rp_a' temp_name '.txt']; % fMRI model specification: Multiple regressors - cfg_files
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
%             matlabbatch{1}.spm.temporal.st.so = [1	7	13	19	25	31	2	8	14	20	26	32	3	9	15	21	27	33	4	10	16	22	28	34	5	11	17	23	29	35	6	12	18	24	30];
        interleaveNum = ceil(sqrt(matlabbatch{k}.spm.temporal.st.nslices));
        nslices = matlabbatch{k}.spm.temporal.st.nslices;
        sliceOrder = [];
        for i=1:interleaveNum
            sliceOrder = [sliceOrder i:interleaveNum:nslices];
        end
         matlabbatch{1}.spm.temporal.st.so = sliceOrder;
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