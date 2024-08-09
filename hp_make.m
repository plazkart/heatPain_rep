%This script is devoted to process all the data obtained with thermal pain
%stimuli

% AYakovlev 13-10-2023
% AYakovlev 13-11-2023

%% Main steps
% mainStruct = hp_make('load'); 

% id= 3;
% mainStruct = hp_make('sp_proc', id, 'sham', 1);
% mainStruct = hp_make('sp_proc', id, 'act', 1);
% hp_make('copyData', id,  'F:\heatPainSP', 'tp_mrs');
% hp_make('copyData', id,  'F:\heatPainSP', 'water_mrs');
% %%

% %%
% mainStruct = hp_make('save', mainStruct);
% %%
% mainStruct =  hp_make('new_field','mainStruct.(nam).proc_check.timepoints_spectra.sham', 0);
% mainStruct =  hp_make('new_field','mainStruct.(nam).proc_check.timepoints_spectra.act', 0);
% %%
% hp_make('getTableData', 'F:\heatPainSP\results', 5)
% %%
% %  mainStruct = hp_make('new', mainStruct);
% %  mainStruct = hp_make('parse', mainStruct);
% %    mainStruct = hp_make('parseMRI', mainStruct);
% %  mainStruct = hp_make('parseFuncTable', mainStruct);
% 
% %  hp_make('save', mainStruct);
%  hp_make('mrs_task_file', 3);
%    hp_make('fmri_proc', 8, mainStruct, 1);
   
%%
% hp_make('expDynamicsSNR'); 
% Make main structure
% if isenv(varargin)
%     hp_make(action, varargin);
% else
%     hp_make(action);
% end
% hp_make(action, varargin);

function [mainStruct, varargout] = hp_make(action, varargin)
% Get into the function varargin as a cell of a number of input values
switch action
    case 'default'
        mainStruct = struct();
        
        mainStruct.meta.date = '13102023';
        mainStruct.meta.subNumbers = 0;
        mainStruct.meta.folder = 'G:\_other\fMRS-heatPain';
        
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
        mainStruct.(nam).proc.dummy_time = 12;
        mainStruct.(nam).proc.slice_order = 'interleaved';
        
        mainStruct.(nam).proc_check.timepoints_spectra.sham = 0;
        mainStruct.(nam).proc_check.timepoints_spectra.act = 0;
        mainStruct.(nam).proc_check.tp_spectra_res.sham = 0;
        mainStruct.(nam).proc_check.tp_spectra_res.act = 0;
        mainStruct.(nam).proc_check.all.act = 0;
        mainStruct.(nam).proc_check.all.sham = 0;
        mainStruct.(nam).proc_check.all.mrs_QA = 1;
        mainStruct.(nam).proc_check.sp_segmented = 0;
        
        mainStruct.meta.subNumbers = mainStruct.meta.subNumbers+1;

        mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\meta']);
        txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'w');
        fprintf(txt_protocol, 'Data processing steps log-file. \n Subject: %s \n Date: %s \n', nam, datetime("today"));
        fclose(txt_protocol);
        
        %% 
        % parsing data and place it into correct directory
    case 'new_field'
        %use as hp_make('new_field', fieldname, default_value);
        % It is needed when, for some reason, there is a need for a new
        % field in te data structure
        mainStruct = hp_make('load'); 
        fieldname = varargin{1};
        default_value = varargin{2};
        for i=1:mainStruct.meta.subNumbers
            nam = sprintf('sub_%02i', i); 
            if ~ischar(default_value)
                eval([fieldname ' = ' num2str(default_value) ';']);
            else
                eval([fieldname ' = default_value;']);
            end

        end
        mainStruct = hp_make('save', mainStruct);

    case 'parse'
        if length(varargin)>0
            mainStruct = varargin{1};
            if length(varargin)>1
                who_id = varargin{2};
            else
                who_id = mainStruct.meta.subNumbers;
            end
        end
        nam = sprintf('sub_%02i', who_id);

        if mainStruct.(nam).data_check.sp<1
            fils_sp = dir([mainStruct.meta.folder mainStruct.(nam).folder '\*.SDAT']);
            if ~isempty(fils_sp)
                mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\sp']);
                for i=1:length(fils_sp)
                    temp = regexp( fils_sp(i).name ,'[_]\d*[_]','match'); %take series number
                    a = strfind(temp{1}, '_');
                    temp{1}(a) = '';
                    seriesNum(i) = str2num(temp{1});                    
                end
                [~, I] = sort(seriesNum);
                caseI = '';

                for i=1:length(I)
                    x = i;
                    switch x
                        case 1
                            mainStruct.(nam).sp.sham = fils_sp(I(i)).name;
                            mainStruct.(nam).data_check.sp = mainStruct.(nam).data_check.sp + 1;
                            caseI = 'sham';
                        case 2
                            if length(I)>3
                                mainStruct.(nam).sp.ref = fils_sp(I(i)).name;
                                mainStruct.(nam).data_check.sp= mainStruct.(nam).data_check.sp + 10;
                                caseI = 'ref';
                            else
                                mainStruct.(nam).sp.act = fils_sp(I(i)).name;
                                mainStruct.(nam).data_check.sp= mainStruct.(nam).data_check.sp + 100;
                                caseI = 'act';
                            end
                        case 4
                            mainStruct.(nam).sp.act = fils_sp(I(i)).name;
                            mainStruct.(nam).data_check.sp= mainStruct.(nam).data_check.sp + 100;
                            caseI = 'act';
                    end
                    copyfile([mainStruct.meta.folder mainStruct.(nam).folder '\' mainStruct.(nam).sp.(caseI)],...
                        [mainStruct.meta.folder mainStruct.(nam).folder '\sp\' nam '_' caseI '.SDAT']);
                    copyfile([mainStruct.meta.folder mainStruct.(nam).folder '\' mainStruct.(nam).sp.(caseI)(1:end-4) 'SPAR'],...
                        [mainStruct.meta.folder mainStruct.(nam).folder '\sp\' nam '_' caseI '.SPAR']);
                end                    
            end
        end
    case 'parseMRI' 
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
                    copyfile([mainStruct.meta.folder mainStruct.(nam).folder '\' fils_mri(i).name],...
                        [mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_func.nii']);
                    mainStruct.(nam).data_check.fmri = 1;
                end
                modalityI = regexp( fils_mri(i).name ,'T2','match');
                if ~isempty(modalityI)
                    copyfile([mainStruct.meta.folder mainStruct.(nam).folder '\' fils_mri(i).name],...
                        [mainStruct.meta.folder mainStruct.(nam).folder '\anat\' nam '_t2check.nii']);
                    mainStruct.(nam).data_check.t2check = 1;
                end
                
            end
        end

    case 'parseFuncTable'
        %use as mainStruct = hp_make('parseFuncTable', mainStruct, id);
        mainStruct = hp_make('load');
        if length(varargin)>0
%             mainStruct = varargin{1};
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
            copyfile([fils_xl(2).folder '\' fils_xl(2).name], [mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
            if length(fils_xl)>2
                mainStruct.(nam).funcTable.fmrs = fils_xl(3).name;
                mainStruct.(nam).data_check.funcTable =mainStruct.(nam).data_check.funcTable+100;
                copyfile([fils_xl(3).folder '\' fils_xl(3).name], [mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_fmrs.xlsx']);
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
        hp_make('save', mainStruct);

    case 'PsychoPy csv file'
        % hp_make('PsychoPy csv file', id, MRcase);
        mainStruct = hp_make('load');
        who_id = varargin{1};
        nam = sprintf('sub_%02i', who_id);
        MRcase = varargin{2};

        fils_csv = dir([mainStruct.meta.folder mainStruct.(nam).folder '\*.csv']);
        for i = 1:length(fils_csv)
            csvVal(i, 1) = fils_csv(i).datenum;
        end
        
        switch MRcase
            case 'MRS'
                [~, lateCSV] = max(csvVal);
                if lateCSV<2
                    nam
                    warning("May be this file is not for MRS");
                end

            case 'MRI'
                [~, lateCSV] = min(csvVal); 
        end
        csv_file = readmatrix([fils_csv(lateCSV).folder '\' fils_csv(lateCSV).name]);
                % get info from csv for MRS
                start_of_MR = csv_file(:, 15);
                start_of_MR = start_of_MR(~isnan(start_of_MR));
                triggerTimes = csv_file(:, 17);
                triggerTimes = triggerTimes(~isnan(triggerTimes));
        
                mainStruct.(nam).proc.(MRcase).startMR = start_of_MR;
                mainStruct.(nam).proc.(MRcase).TTLtimes = triggerTimes;
        
                %get pain estimation from the file
                try
                 [rating, reaction_time] = AssessmentParser2([fils_csv(lateCSV).folder '\' fils_csv(lateCSV).name]);
                catch ME
                    if contains(ME.message, 'key_resp_2_rt')
                        return;
                        rating = 0;
                        reaction_time = 0;
                    end
                end
                       
                mainStruct.(nam).proc.(MRcase).rating = rating;
                mainStruct.(nam).proc.(MRcase).reaction_time = reaction_time;

        hp_make('save', mainStruct);


%% input-output tasks

    case 'save'
        if length(varargin)>0
            mainStruct = varargin{1};
        end
        mainStruct.meta.date = datetime("today");
        save([mainStruct.meta.folder '\_meta\structMeta.mat'], 'mainStruct');
        
        
    case 'load'
        hp_default;
%         load([mainStruct.meta.folder '\_meta\structMeta.mat']);
%% Processing steps
    %fmri
     
    
    case 'fmri_proc'
        % use it as hp_make('fmri_proc', id, mainStruct, procSteps); mainStruct - if
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

        if mainStruct.(nam).data_check.fmri == 1
            fmri_img = spm_vol([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_func.nii']);
            NSA = length(fmri_img);
            for i=1:NSA
                func_data(i, 1) = {[mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_func.nii,' num2str(i)]};
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
            callfMRIProcessing(inputs, mainStruct.(nam).proc.slice_order, 1);


            fils_func = dir([mainStruct.meta.folder mainStruct.(nam).folder '\func\*' nam '*']);
            mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\derived']);
            for i=1:length(fils_func)
                if ~contains(fils_func(i).name, 'xlsx')
                    if ~strcmp(fils_func(i).name, [nam '_func.nii'])
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
            for i=1:NSA
                func_data(i, 1) = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\swra' nam '_func.nii,' num2str(i)]};
            end
            inputs{2, crun} = func_data; % fMRI model specification: Scans - cfg_files
            if floor(mod(mainStruct.(nam).data_check.funcTable, 100)/10)>0
                regressorList = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
%                 regressorList(:, 1) = regressorList(:, 1) - regressorList(1, 1) - (12 - mainStruct.(nam).proc.dummy_time);
                TTL_time = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
                regressorList(:, 1) = regressorList(:, 1) - TTL_time(1) - (12 - mainStruct.(nam).proc.dummy_time);
            end
            inputs{3, crun} = regressorList(:, 1); % fMRI model specification: Onsets - cfg_entry
            inputs{4, crun} = regressorList(:, 2); % fMRI model specification: Durations - cfg_entry
            inputs{5, crun} = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\rp_a' nam '_func.txt']}; % fMRI model specification: Multiple regressors - cfg_files
        end
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs, inputs{:});

        % MRI processing
    case 'T2_movements_check'
        % use as hp_make('T2_movements_check', id)
        mainStruct = hp_make('load');

        if nargin<2
            error('There no subject name to process');
        end
        id = varargin{1};
        nam = sprintf('sub_%02i', id);

        %look for t2 images
        filst2 = dir([mainStruct.meta.folder mainStruct.(nam).folder '\*T2*.nii']);
        if length(filst2)>1
            for i =1:2
%                 series = split(filst2(i).name, '_');
%                 series = regexp(series{end}, '\d*', 'match');
%                 series = series{end};
%                 series = str2num(series);
%                 filst2(i).series = series;
                input2align{i, 1} = [filst2(i).folder '\' filst2(i).name];
            end

            alignImages(input2align);
            filtxt = dir([mainStruct.meta.folder mainStruct.(nam).folder '\rp_*.txt']);
            copyfile([filtxt.folder '\' filtxt.name], [filtxt.folder '\meta\displacement_t2.txt']);
            filtxt = fopen([filtxt.folder '\meta\displacement_t2.txt'], 'r');
            A = fscanf(filtxt, '%f');
            A = A(7:9);
            A = sqrt(sum(A.^2));
            mainStruct.(nam).proc.displacentT2 = A;

            mainStruct = hp_make('save', mainStruct);

        end

            
  

% MRS processing
        
    case 'sp_proc'
        %% use as hp_make('sp_proc', id, sp_name, time_points_action)
        % id - id of the precossed subject
        % sp_name - type of the processed spectra {'sham', 'ref', 'act'}
        % time_points_action - is it needed to process as time points? (0
        % no, 1 - yes)
        mainStruct = hp_make('load');

        if nargin<2
            error('There no subject name to process');
        end
        id = varargin{1};
        nam = sprintf('sub_%02i', id);
        
        cases = {'sham', 'ref', 'act'};
        sp_name = varargin{2};
        sp_action = varargin{3};
        k=0;
        for i=1:length(cases)
            if strcmp(sp_name, cases{i})
                k=i;
            end
        end
        if k==0
            error('there is no such spectral state (only sham, ref and act)');
        else
            if ~floor(mod(mainStruct.(nam).data_check.funcTable, 10^k)/10^(k-1))>0 %this takes a k-th number in data check field
                error('there is no such data');
            end
        end
        
       switch sp_action
           case 0
               [mainStruct, sp_out] = hp_make('spectraPreprocessing', [mainStruct.meta.folder '\' nam '\sp\' nam '_' sp_name '.SDAT'], mainStruct, id, sp_name);
               sp_fpa_av_wr = sp_out;

               % process as whole
               txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'a');
               new_dir = [mainStruct.meta.folder '\' nam '\sp\derived\'];
               mkdir(new_dir);
               sp_naming = sprintf('%s_all_%s', nam, sp_name);
               mrs_writeSDAT([new_dir sp_naming '.SDAT'],...
                   sp_fpa_av_wr.fids);
               copyfile([mainStruct.meta.folder '\' nam '\sp\' nam '_' sp_name '.SPAR'],...
                   [new_dir sp_naming '.SPAR']);
               my_editSPAR([new_dir sp_naming '.SPAR'],{'rows','spec_num_row', 'spec_row_upper_val'}, {1,1,1});
               fprintf( txt_protocol, 'saved in: %s \n', [new_dir sp_naming '.SDAT']);
               fclose(txt_protocol);

           case 1
               %process as interleaved time points
               %be sure that you have a time-points table
               if isfield(mainStruct.(nam).proc, 'tp_matrix')
                   k = max(mainStruct.(nam).proc.tp_matrix); % maximal index of the time point
               else
                   error("There is no time points matrix for data grouping. Make 'mrs_task_file' case");
               end
               sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\' nam '_' sp_name '.SDAT'], 1);
               for i=1:k
                   time_point_dynamics = [1:315];
                   time_point_dynamics = time_point_dynamics(mainStruct.(nam).proc.tp_matrix==i);
                   time_point_dynamics = time_point_dynamics(time_point_dynamics<sp.sz(2));
                   sp_tp = op_takeaverages(sp, time_point_dynamics);
                   tp_nam = sprintf('tp_%02i', i);
                   [mainStruct, sp_out] = hp_make('spectraPreprocessing', sp_tp , mainStruct, id, [tp_nam '_' sp_name]);
                   %                 [mainStruct, sp_out] = hp_make('spectraPreprocessing', ...
                   %                 sp_tp , mainStruct, id, [tp_nam '_' sp_name], [4.2 5 1]);
                   %                 %Made once for sub_11 /19-01-2024
                   sp_fpa_av_wr = sp_out;

                   % save into the protocol
                   txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'a');
                   new_dir = [mainStruct.meta.folder '\' nam '\sp\derived\'];
                   mkdir(new_dir);
                   sp_naming = sprintf('%s_all_%s', nam, [tp_nam '_' sp_name]);
                   mrs_writeSDAT([new_dir sp_naming '.SDAT'],...
                       sp_fpa_av_wr.fids);
                   copyfile([mainStruct.meta.folder '\' nam '\sp\' nam '_' sp_name '.SPAR'],...
                       [new_dir sp_naming '.SPAR']);
                   my_editSPAR([new_dir sp_naming '.SPAR'],{'rows','spec_num_row', 'spec_row_upper_val', 'samples', 'spec_num_col'}, {1,1,1,1024*16, 1024*16});
                   fprintf( txt_protocol, 'saved in: %s \n', [new_dir sp_naming '.SDAT']);
                   fclose(txt_protocol);
               end
               mainStruct.(nam).proc_check.timepoints_spectra.(sp_name) = 1;
               mainStruct = hp_make('save', mainStruct);
           case 2
               sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\' nam '_' sp_name '.SDAT'], 1);
               for i=1:3
                   if 100*i>sp.averages
                       sp_p{i} = op_takeaverages(sp, 100*i-99:sp.averages);
                   else
                       sp_p{i} = op_takeaverages(sp, 100*i-99:100*i);
                   end
                   [~, sp_out] = hp_make('spectraPreprocessing', sp_p{i}, mainStruct, id, sp_name);
                   sp_fpa_av_wr = sp_out;
                   sp_naming = sprintf('%s_p%i_%s', nam, i, sp_name);
                    new_dir = [mainStruct.meta.folder '\' nam '\sp\derived\'];
                   mrs_writeSDAT([new_dir sp_naming '.SDAT'],...
                       sp_fpa_av_wr.fids);
                   copyfile([mainStruct.meta.folder '\' nam '\sp\' nam '_' sp_name '.SPAR'],...
                       [new_dir sp_naming '.SPAR']);
                   my_editSPAR([new_dir sp_naming '.SPAR'],{'rows','spec_num_row', 'spec_row_upper_val'}, {1,1,1});
                   my_editSPAR([new_dir sp_naming '.SPAR'],{'samples'}, {2048*16});
                   txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'a');
                   fprintf( txt_protocol, 'saved in: %s \n', [new_dir sp_naming '.SDAT']);
                   fclose(txt_protocol);

               end

           case 3
               %process as interleaved time points
               %added time smoothing
               if isfield(mainStruct.(nam).proc, 'tp_matrix')
                   k = max(mainStruct.(nam).proc.tp_matrix); % maximal index of the time point
               else
                   error("There is no time points matrix for data grouping. Make 'mrs_task_file' case");
               end
               sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\' nam '_' sp_name '.SDAT'], 1);
               %make smoothin data
               sp.fids = 0.5*sp.fids + 0.25*circshift(sp.fids, 1, 1)+0.25*circshift(sp.fids, -1, 1);
               sp.specs = 0.5*sp.specs + 0.25*circshift(sp.specs, 1, 1)+0.25*circshift(sp.specs, -1, 1);
               for i=1:k
                   time_point_dynamics = [1:315];
                   time_point_dynamics = time_point_dynamics(mainStruct.(nam).proc.tp_matrix==i);
                   time_point_dynamics = time_point_dynamics(time_point_dynamics<sp.sz(2));
                   sp_tp = op_takeaverages(sp, time_point_dynamics);
                   tp_nam = sprintf('tp_%02i', i);
                   [mainStruct, sp_out] = hp_make('spectraPreprocessing', sp_tp , mainStruct, id, [tp_nam '_' sp_name]);
                   %                 [mainStruct, sp_out] = hp_make('spectraPreprocessing', ...
                   %                 sp_tp , mainStruct, id, [tp_nam '_' sp_name], [4.2 5 1]);
                   %                 %Made once for sub_11 /19-01-2024
                   sp_fpa_av_wr = sp_out;

                   % save into the protocol
                   txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'a');
                   new_dir = [mainStruct.meta.folder '\' nam '\sp\derived\'];
                   mkdir(new_dir);
                   sp_naming = sprintf('%s_%s_sm', nam, [tp_nam '_' sp_name]);
                   mrs_writeSDAT([new_dir sp_naming '.SDAT'],...
                       sp_fpa_av_wr.fids);
                   copyfile([mainStruct.meta.folder '\' nam '\sp\' nam '_' sp_name '.SPAR'],...
                       [new_dir sp_naming '.SPAR']);
                   my_editSPAR([new_dir sp_naming '.SPAR'],{'rows','spec_num_row', 'spec_row_upper_val', 'samples', 'spec_num_col'}, {1,1,1,1024*16, 1024*16});
                   fprintf( txt_protocol, 'saved in: %s \n', [new_dir sp_naming '.SDAT']);
                   fclose(txt_protocol);
               end
               mainStruct.(nam).proc_check.timepoints_spectra.smoothed = 1;
               mainStruct = hp_make('save', mainStruct);

       


           
        end
    
    case 'LinewidthAssessment'
        %use as hp_make('LinewidthAssessment', id, condition, tp, bc)
        % assess FWHM of Cr and NAA signals
        mainStruct = hp_make('load');
        id = varargin{1};
        condition = varargin{2};
        tp = varargin{3};
        if length(varargin)>3
            bc = varargin{4};
        else 
            bc = 0;
        end

        nam = sprintf('sub_%02i', id);
        
        switch tp
            case 1
                for i=1:6
                    sp_nam = sprintf('tp_%02i', i);
                    if bc
                        sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\derived\' nam '_all_' sp_nam '_' condition '_bc.SDAT'], 1);
                    else
                        sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\derived\' nam '_all_' sp_nam '_' condition '.SDAT'], 1);
                    end
                    %% assess linewidth
                    Estimates = MRSassessment(sp);
                    mainStruct.(nam).proc.(condition).(sp_nam) = Estimates;
                    if bc
                        mainStruct.(nam).proc.(condition).(sp_nam).LWCr_bc = Cr_fwhm;
                        mainStruct.(nam).proc.(condition).(sp_nam).LWNAA_bc = NAA_fwhm;
                    else
%                         mainStruct.(nam).proc.(condition).(sp_nam).LWCr = Cr_fwhm;
%                         mainStruct.(nam).proc.(condition).(sp_nam).LWNAA = NAA_fwhm;
%                         mainStruct.(nam).proc.(condition).(sp_nam).HCr = Cr_H;
%                         mainStruct.(nam).proc.(condition).(sp_nam).HNAA = NAA_H;
                    end

                end
            case 3
                %temporally smoothed_data
                for i=1:6
                    sp_nam = sprintf('tp_%02i', i);
                    if bc
                        sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\derived\' nam '_all_' sp_nam '_' condition '_bc.SDAT'], 1);
                    else
                        sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\derived\' nam '_' sp_nam '_' condition '_sm.SDAT'], 1);
                    end
                    %% assess linewidth
                    Estimates = MRSassessment(sp);
                    sp_nam = sprintf('tp_%02i_sm', i);
                    mainStruct.(nam).proc.(condition).(sp_nam) = Estimates;
                end
            case 0
                if mainStruct.(nam).proc_check.all.(condition)<1
                    return;
                end
                sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\derived\' nam '_all_' condition '.SDAT'], 1);
                sp.flags.averaged = 1;
                sp.dims.averages = 0;
                sp = op_autophase(sp, 1.8, 3.5, 0);
                Cr_fwhm = op_getLW(sp, 2.8, 3.1);
                NAA_fwhm = op_getLW(sp, 1.8, 2.15);
                mainStruct.(nam).proc.(condition).all.LWCr = Cr_fwhm;
                mainStruct.(nam).proc.(condition).all.LWNAA = NAA_fwhm;
            case 2
                for i=1:3
                    tp_nam = sprintf('p%i', i);
                    sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\derived\' nam '_' tp_nam '_' condition '.SDAT'], 1);
                    sp.flags.averaged = 1;
                    sp.dims.averages = 0;
                    sp = op_autophase(sp, 1.8, 3.5, 0);
                    Cr_fwhm = op_getLW(sp, 2.8, 3.1);
                    NAA_fwhm = op_getLW(sp, 1.8, 2.15);
                    mainStruct.(nam).proc.(condition).(tp_nam).LWCr = Cr_fwhm;
                    mainStruct.(nam).proc.(condition).(tp_nam).LWNAA = NAA_fwhm;
                end
        end
        hp_make('save', mainStruct);

    case 'spectraPreprocessing'
        %use as hp_make('spectraPreprocessing', path, mainStruct, id)
        if ~isstruct(varargin{1})
            path = varargin{1}; 
            sp = io_loadspec_sdat(path, 1);
        else
            sp = varargin{1};
        end                                                                     
        mainStruct = varargin{2};
        id = varargin{3};
        sp_name = varargin{4};
        if length(varargin)<5
            pars{1, 1} = [1.9 2.1 1];
            pars{2, 1} = {[4.1 5.1], 6};
        else
            pars{1, 1} = varargin{5};
        end


        
        nam = sprintf('sub_%02i', id);
        % Processing steps
        % test 1: Freq and Phase aligning using ALL spectra (1) and divided
        % series (2)
        % frequency and phase aligning (NAA signal as ref)
        txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'a');
        fprintf(txt_protocol, '------ \n spectra processing steps: case - %s \n', sp_name);
        fprintf(txt_protocol, '%s \n', datetime("today"));
        pars{3, 1} = 16;
        sp_zf = op_zeropad(sp, pars{3, 1});
        fprintf(txt_protocol, 'op_zeropad: pars - %02i \n', pars{3, 1}(:));

        pars{1, 1} = [1.9 2.1 1]; 
        sp_fpa_zf = op_freqAlignAverages_fd(sp_zf, pars{1, 1}(1), pars{1, 1}(2), pars{1, 1}(3), 'n');
        fprintf(txt_protocol, 'op_freqAlignAverages_fd: pars - %3.1f %3.1f %3.1f \n', pars{1, 1}(:));
        
        sp_fpa_av = op_averaging(sp_fpa_zf );  
        fprintf(txt_protocol, 'op_averaging \n');
        if size(pars, 1)>1
            sp_fpa_av_wr = op_removeWater(sp_fpa_av, pars{2, 1}{1}, pars{2, 1}{2});
            fprintf( txt_protocol, 'op_removeWater: pars - %3.1f %3.1f %02i \n', pars{2, 1}{:});
        else
            sp_fpa_av_wr = sp_fpa_av;
        end
        fclose(txt_protocol);
        varargout{1} = sp_fpa_av_wr;

    case 'mrs_task_file'
        %use as hp_make('mrs_task_file', id)
        mainStruct = hp_make('load');
        %using MEDOC TSA data
        if nargin<2
            error('There no subject name to process');
        end 
        id = varargin{1};
        nam = sprintf('sub_%02i', id);
        k=2;

        %determination of main values
        dummy_Delta = mainStruct.(nam).proc.TTLtimes(1) - mainStruct.(nam).proc.startMR; %time between 1st TSA trigger and start of MR-protocol

        if floor(mod(mainStruct.(nam).data_check.funcTable, 10^k)/10^(k-1))>0
            ttlTime = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_fmrs.xlsx']);
            regressorList = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_fmrs.xlsx']);
            task_starts = regressorList(:,1)-ttlTime(1)+dummy_Delta; % times of 90% of stimulus temperature
            %relative to the 1st dynamic of MRS
            stim_dur = regressorList(:,2); %time between 90% of temperature (star and end)
        else 
            error("There is no stimulation table data to make time_point matrix");
        end

        sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\' nam '_act.SDAT'], 1);
        mrs_NSAmax = sp.averages;
        mrs_timings = [1:mrs_NSAmax]*2-2;
        time_point_matrix = zeros(mrs_NSAmax, 1);
        shifts = zeros(length(task_starts));

        for i=1:length(task_starts)
            tempMRSTimings = mrs_timings - task_starts(i);
            [shifts(i), idx] = min(abs(tempMRSTimings));
            time_point_matrix(idx) = 1;
        end       
        %need to write first active time point
        task_starts(task_starts>mrs_NSAmax*2)=[];
        for k=2:6
            tps = find(time_point_matrix==k-1)+1;
            time_point_matrix(setdiff(tps, find(time_point_matrix==1))) = k;
        end
        time_point_matrix(mrs_NSAmax+1:end) = [];

        %test of stimulus fuction
        if true
            xlFile = readtable([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_fmrs.xlsx']);
            
            
            temp = xlFile.Tec_C_;
            time = xlFile.Timestamp_msec_/1000;
            stimTemp = median(temp(1:1000))+ (max(temp) - median(temp(1:1000)))*0.9;

            time = time-ttlTime(1)+dummy_Delta;
            
            figure(); plot(time, temp); hold on;
            for k=1:6
                tps = find(time_point_matrix==k);
                scatter(mrs_timings(tps), ones(size(tps))*stimTemp);
            end
            saveas(gcf,[mainStruct.meta.folder '\' nam '\meta\MRStiming.jpg'],'jpg'); hold off;
            

        end
        
        mainStruct.(nam).proc.start_dynamic = dummy_Delta;
        mainStruct.(nam).proc.tp_matrix = time_point_matrix;
        mainStruct.(nam).proc.shifts = shifts;
        mainStruct.(nam).proc_check.tp_matrix = 1;

        mainStruct = hp_make('save', mainStruct);
        
    case 'Make_HRFforMRSdata'
        %use as hp_make('Make_HRFforMRSdata', id)
        %make personal HRF for heat pain stimuli during fMRS
        %simple HRF convolved with stimulus function (var task_starts)
        %get standard HRF

        mainStruct = hp_make('load');
        id = varargin{1};
        nam = sprintf('sub_%02i', id);

        time_point_matrix = mainStruct.(nam).proc.tp_matrix;
        
        dummyTime = mainStruct.(nam).proc.TTLtimes(1) - mainStruct.(nam).proc.startMR;
        ttlTime = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_fmrs.xlsx']);
            regressorList = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_fmrs.xlsx']);
            task_starts = regressorList(:,1)-ttlTime(1);  
            stim_dur = regressorList(:,2);
        task_starts = task_starts + dummyTime;

        MR = 16;
        dt = [0:315*MR*2]/MR;
        u = zeros(length(dt), 1);

        %trying to find timestamps in which stimulation was ON
        for i=1:length(task_starts)
%             if i==59
%                 test = 1;
%             end
            a = find(dt>task_starts(i));
            if isempty(a)
                continue
            end
            lower_side = a(1)+ceil(dt(a(1))-task_starts(i));
            a = find(dt>task_starts(i)+3); %Fix stimulus
            if isempty(a)
                continue
            end
            upper_side = a(1)+ceil(dt(a(1))-(task_starts(i)+3));
            u(lower_side:upper_side, 1)=1;
        end

        pars = [6, 16, 1, 1, 6, 0, 32];
        [bf, p] = spm_hrf(1/16, pars(1:6), 16);

        hrf = conv(u, bf);

        %test of hrf
        if true
             hrf = hrf(1:length(dt)); quantileHRF_90 = quantile(hrf, 0.9);
            mrs_timings = [1:315]*2-2;
            figure(); plot(dt, hrf); hold on;
            for k=1:6
                tps = find(time_point_matrix==k);
                scatter(mrs_timings(tps), ones(size(tps))*quantileHRF_90);
            end
            saveas(gcf,[mainStruct.meta.folder '\' nam '\meta\MRS_activation.jpg'],'jpg'); hold off;
            

        end

        TR_samples = [0:315]*2*MR+1;
        hrf_samples=hrf(TR_samples);
        hrf_samples(1)=[];
        %% temporal smoothing
         hrf_smoothed = 0.5*hrf+0.25*circshift(hrf, 2*MR)+0.25*circshift(hrf, -2*MR);
         hrf_samples_smoothed=hrf_smoothed(TR_samples);hrf_samples_smoothed(1)=[];
        for i=1:max(time_point_matrix)
            hrf_mean(i) = mean(hrf_samples(time_point_matrix==i));
            hrf_mean_smoothed(i) = mean(hrf_samples_smoothed(time_point_matrix==i));
        end

        mainStruct.(nam).proc.hfr_mrs = hrf_mean;
        mainStruct.(nam).proc.hfr_mrs_smoothed = hrf_mean_smoothed;
        mainStruct = hp_make('save', mainStruct);

    case 'copyData'
        %use as hp_make('copyData', id,  toDir, type)
        %where id  - is num of subject, toDir - where to copy files, type -
        %what kind of data should be copied (see below in the switch cases)
        id = varargin{1};
        nam = sprintf('sub_%02i', id);
        out_path = varargin{2};
        type = varargin{3};
        mainStruct = hp_make('load');
        switch type
            case 'all_mrs'
                cases = {'sham', 'act'};
                for ii =1:2
                    copyFiles = dir([mainStruct.meta.folder '\' nam '\sp\derived\*all_' cases{ii} '.*']);
                    for i = 1:length(copyFiles)
                        if ~contains(copyFiles(i).name, 'bc')
                            copyfile([copyFiles(i).folder '\' copyFiles(i).name], [out_path '\' copyFiles(i).name]);
                        end
                    end
                end
                hp_make('copyData', id,  out_path, 'makeProcessingList')

            case 'sixPoints_mrs'
                cases = {'sham', 'act'};
                for ii =1:2
                    copyFiles = dir([mainStruct.meta.folder '\' nam '\sp\derived\*_p*' cases{ii} '.*']);
                    for i = 1:length(copyFiles)
                        if ~contains(copyFiles(i).name, 'bc')
                            copyfile([copyFiles(i).folder '\' copyFiles(i).name], [out_path '\' copyFiles(i).name]);
                        end
                    end
                end
                hp_make('copyData', id,  out_path, 'makeProcessingList')

            case 'tp_mrs'
                if mainStruct.(nam).proc_check.timepoints_spectra.sham<1
                    error('There is no time points divided data YET');
                end
                %in case of initial data
                copyFiles = dir([mainStruct.meta.folder '\' nam '\sp\derived\*tp*']);
                %in case of smoothed data
                copyFiles = dir([mainStruct.meta.folder '\' nam '\sp\derived\*tp*sm*']);
                for i = 1:length(copyFiles)
                    if ~contains(copyFiles(i).name, 'bc')
                        copyfile([copyFiles(i).folder '\' copyFiles(i).name], [out_path '\' copyFiles(i).name]);
                    end
                end
                hp_make('copyData', id,  out_path, 'makeProcessingList')
            case 'tp_mrs_bc'
                if mainStruct.(nam).proc_check.timepoints_spectra.sham<1
                    error('There is no time points divided data YET');
                end
                if mainStruct.(nam).proc_check.bold_correction<1
                    error('There is no BOLD correction done yet');
                end
                copyFiles = dir([mainStruct.meta.folder '\' nam '\sp\derived\*tp*bc*']);
                for i = 1:length(copyFiles)
                    copyfile([copyFiles(i).folder '\' copyFiles(i).name], [out_path '\' copyFiles(i).name]);
                end
                hp_make('copyData', id,  out_path, 'makeProcessingList')

            case 'water_mrs'
                k = 2;
                if ~floor(mod(mainStruct.(nam).data_check.sp, 10^k)/10^(k-1))>0 %this takes a k-th number in data check field
                    error('there is no water reference SDAT-file or no info in mainStruct');
                end
                copyFiles = dir([mainStruct.meta.folder mainStruct.(nam).folder '\sp\' nam '_ref.*']);
                for i = 1:length(copyFiles)
                    copyfile([copyFiles(i).folder '\' copyFiles(i).name], [out_path '\' copyFiles(i).name]);
                end

            case 'makeProcessingList'
                fils = dir([ out_path '\*' nam '*.SDAT' ]);
                fils_ref = dir([ out_path '\*' nam '_ref.SDAT' ]);
%                 fils_lists = dir([out_path '\listProc.txt' ]);
                if isempty(fils)
                    error('There is no files to make processing list');
                end
                txtFil = fopen([out_path '\listProc.txt' ], 'a');

                for i = 1:length(fils)
                    if ~contains(fils(i).name, 'ref')
                        fprintf(txtFil, '%s %s\n', fils(i).name, fils_ref(1).name);
                    end
                end
                fclose(txtFil);

            





        end

    case 'estimations_parse'
        %use as hp_make('estimations_parse', id)
        hp_default
        id = varargin{1};
        nam = sprintf('sub_%02i', id);
        if ~isempty(mainStruct.(nam).funcTable.est)
            [rating, reaction_time] = AssessmentParser2([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_est.csv']);
            mainStruct.(nam).proc.react_time = reaction_time;
            mainStruct.(nam).proc.rating = rating;
        end
        mainStruct = hp_make('save', mainStruct);
            

        %% Results parsing and summarising
    case 'getTableData2'
        %use as hp_make('getTableData2', path, id, tp_case)
        %Since there is more cases than before, a new release of the
        %function is needed
        %tp_case: if 2 - for six points processing pipeline
        mainStruct = hp_make('load');
        res_path = varargin{1};
        id = varargin{2};
        nam = sprintf('sub_%02i', id);
        fils = dir([res_path '\' nam '*\table']);
        if length(varargin)>2
            tp_case = varargin{3};
        else
            tp_case = 0;
        end
        
        switch tp_case
            case 2
                for i=1:length(fils)
                    temp = split(fils(i).folder, '_');
                    temp{end-1} = regexp(temp{end-1}, '\d', 'match');
                    tp_num = str2num(temp{end-1}{1});
                    mod_case = temp{end};

                    out_dir = sprintf('p%i_%s',tp_num, mod_case);
                    mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir]);
                    copyfile([fils(i).folder '\' fils(i).name], [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\' fils(i).name]);
                    copyfile([fils(i).folder '\coord'], [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\coord']);
                    
                    tp_nam = sprintf('p%i', tp_num);
                    mainStruct.(nam).proc.(mod_case).(tp_nam).exist = 1;
                    mainStruct.(nam).proc.(mod_case).(tp_nam).path = [mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir '\' fils(i).name];


                    mainStruct = hp_make('processLCTable',mainStruct, id, mod_case, tp_nam);
                end
                case 3
                    %for smoothed data
                for i=1:length(fils)
                    temp = split(fils(i).folder, '_');
                    tp_num = str2num(temp{end-2});
                    mod_case = temp{end-1};

                    out_dir = sprintf('tp_%02i_%s_sm',tp_num, mod_case);
                    mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir]);
                    copyfile([fils(i).folder '\' fils(i).name], [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\' fils(i).name]);
                    copyfile([fils(i).folder '\coord'], [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\coord']);
                    
                    tp_nam = sprintf('tp_%02i_sm', tp_num);
                    mainStruct.(nam).proc.(mod_case).(tp_nam).exist = 1;
                    mainStruct.(nam).proc.(mod_case).(tp_nam).path = [mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir '\' fils(i).name];


                    mainStruct = hp_make('processLCTable',mainStruct, id, mod_case, tp_nam);
                end
        mainStruct = hp_make('save', mainStruct);


        end

    case 'getTableData'
        %use as hp_make('getTableData', path, id, tp_case)
        mainStruct = hp_make('load');
        res_path = varargin{1};
        id = varargin{2};
        nam = sprintf('sub_%02i', id);
        fils = dir([res_path '\' nam '*\table']);
        if length(varargin)>2
            tp_case = varargin{3};
        else
            tp_case = 0;
        end
        for i=1:length(fils)
            temp = split(fils(i).folder, '_');
            if contains(fils(i).folder, 'bc')
                bc_case = 1;
                temp = temp(1:end-1);
            else 
                bc_case = 0;
            end
            mod_case = temp{end};
            if contains(fils(i).folder, '_tp_')
                tp_case = 1;
                tp_num = str2num(temp{end-1});
            end            
            if tp_case
                if ~bc_case
                    out_dir = sprintf('tp_%02i_%s', tp_num, mod_case);
                else
                    out_dir = sprintf('tp_%02i_%s_bc', tp_num, mod_case);
                end
            else
                out_dir = sprintf('all_%s', mod_case);
            end
            mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir]);
            copyfile([fils(i).folder '\' fils(i).name], [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\' fils(i).name]);
            copyfile([fils(i).folder '\coord'], [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\coord']);
            
            if tp_case
                if ~bc_case
                    tp_nam = sprintf('tp_%02i', tp_num);
                else
                    tp_nam = sprintf('tp_%02i_bc', tp_num);
                end
                mainStruct.(nam).proc.(mod_case).(tp_nam).exist = 1;
                mainStruct.(nam).proc.(mod_case).(tp_nam).path = [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\' fils(i).name];
                mainStruct = hp_make('processLCTable',mainStruct, id, mod_case, tp_nam);
                mainStruct.(nam).proc_check.tp_spectra_res.(mod_case) = tp_num;
            else
                mainStruct.(nam).proc.(mod_case).all.exist = 1;
                mainStruct.(nam).proc.(mod_case).all.path = [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\' fils(i).name];
                mainStruct = hp_make('processLCTable',mainStruct, id, mod_case, 'all');
%                 mainStruct.(nam).proc_check.tp_spectra_res.(mod_case) = tp_num;
            end
        end
         mainStruct = hp_make('save', mainStruct);

    case 'processLCTable'
        %use as hp_make('processLCTable',mainStruct, id, condition, tp_number)
        mainStruct =  varargin{1};
        id =  varargin{2};
        condition =  varargin{3};
        tp_nam =  varargin{4};
        nam = sprintf('sub_%02i', id);
        

        resTab = mrs_readLcmodelTABLE(mainStruct.(nam).proc.(condition).(tp_nam).path);
        mainStruct.(nam).proc.(condition).(tp_nam).NAA = resTab.concentration(20);
        mainStruct.(nam).proc.(condition).(tp_nam).Cr = resTab.concentration(21);
        mainStruct.(nam).proc.(condition).(tp_nam).Glx = resTab.concentration(22);
        mainStruct.(nam).proc.(condition).(tp_nam).NAAerr = resTab.SDpct(20);
        mainStruct.(nam).proc.(condition).(tp_nam).Crerr = resTab.SDpct(21);
        mainStruct.(nam).proc.(condition).(tp_nam).Glxerr = resTab.SDpct(22);
%% MRS quality control procedures
% draw MR-spectra
    case 'coord_drawSpectra'
        %use as hp_make('coord_drawSpectra', id, condition, tp_number)
        mainStruct = hp_make('load');
        condition = varargin{2};
        id = varargin{1};
        nam = sprintf('sub_%02i', id);

        if mainStruct.(nam).proc_check.all.(condition)>0
            coordPath = [mainStruct.meta.folder mainStruct.(nam).folder  '\results\sp\all_' condition '\coord'];
            sp_figure = io_readlcmcoord_getBackground(coordPath,'sp');
            fit_figure = io_readlcmcoord_getBackground(coordPath,'fit');
            plot(sp_figure.ppm, [sp_figure.specs, fit_figure.specs]);
            set(gca,'Xdir','reverse');
            exportgraphics(gcf,[mainStruct.meta.YDfolder '\spectraFigures\' nam '_' condition '.png'], 'BackgroundColor', 'white');
        end


    case 'spVoxelPlacement'
        % use as hp_make('spVoxelPlacement', id, condition)
        %makes voxel mask in MNI-space
        mainStruct = hp_make('load');
        id = varargin{1};
        condition = varargin{2};
        nam = sprintf('sub_%02i', id);

        MRS_struct = CoRegStandAlone({[mainStruct.meta.folder mainStruct.(nam).folder '\sp\' nam '_' condition '.SDAT']},...
            {[mainStruct.meta.folder mainStruct.(nam).folder '\anat\' nam '_anat.nii']});
        mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\anat\derived']);
        copyfile(MRS_struct.mask.vox1.outfile{1}, [mainStruct.meta.folder mainStruct.(nam).folder '\anat\derived\' nam '_sp_mask_' condition '.nii']);
        delete(MRS_struct.mask.vox1.outfile{1});
        
        mainStruct.(nam).proc.(condition).all.path = [mainStruct.meta.folder mainStruct.(nam).folder '\anat\derived\' nam '_sp_mask_' condition '.nii'];
        mainStruct.(nam).proc.(condition).all.GM = MRS_struct.out.vox1.tissue.fGM;
        mainStruct.(nam).proc.(condition).all.WM = MRS_struct.out.vox1.tissue.fWM;
        mainStruct.(nam).proc.(condition).all.CSF = MRS_struct.out.vox1.tissue.fCSF;

        mainStruct.(nam).proc_check.sp_segmented = 1;
        mainStruct = hp_make('save', mainStruct);
        %translate to MNI
        callNormilise(mainStruct.(nam).proc.(condition).all.path,...
            [mainStruct.meta.folder mainStruct.(nam).folder '\anat\y_' nam '_anat.nii']);
        
    case 'count_fmriMetrics'
        % use as [~, metric] = hp_make('count_fmriMetrics', id, metrics)
        % 'beta' - quantifying mean ratio of beta coef. 1 (stim) and
        % constant term (last one)

        %Check this for Visual stimulation data
        mainStruct = hp_make('load');
        id = varargin{1};
        metrics = varargin{2};
        nam = sprintf('sub_%02i', id);

%         [filepath,name] = fileparts(mainStruct.(nam).proc.sp_mask.path);
            sp_mask_path = [mainStruct.meta.folder '\' nam '\anat\derived\w' nam '_sp_mask_sham.nii'];
           callSmothing(sp_mask_path);
         Vox_mask = spm_vol([mainStruct.meta.folder '\' nam '\anat\derived\sw' nam '_sp_mask_sham.nii']);
        
        switch metrics
            case 'beta'
                fils = dir([mainStruct.meta.folder '\' nam '\derived\res\beta_*']);
                img_beta1= [fils(1).folder '\' fils(1).name];
                img_beta2 = [fils(end).folder '\' fils(end).name];
                V_img_1 = spm_vol(img_beta1);
                makeSameResolution(V_img_1.fname, Vox_mask.fname);
                img_mask = [filepath '\rsw' name '.nii'];

                BOLD_beta = countBeta(img_mask, img_beta1, img_beta2);
                mainStruct.(nam).proc.bold.mean_delta_beta1 = BOLD_beta;
                varargout{1} = BOLD_beta;
                hp_make('save', mainStruct);

            case 'insula_beta'
                fils = dir([mainStruct.meta.folder '\' nam '\derived\res\beta_*']);
                img_mask = [mainStruct.meta.folder '\_meta\atlas_map\rsinsula_atlas.nii'];
                img_beta1= [fils(1).folder '\' fils(1).name];
                img_beta2 = [fils(end).folder '\' fils(end).name];
                 BOLD_beta = countBeta(img_mask, img_beta1, img_beta2);
                varargout{1} = BOLD_beta;

            case 'cluster'
                fils = dir([mainStruct.meta.folder '\' nam '\derived\res\beta_*']);
                img_mask = [fils(1).folder '\clus1.nii'];
                img_beta1= [fils(1).folder '\' fils(1).name];
                img_beta2 = [fils(end).folder '\' fils(end).name];
                 BOLD_beta = countBeta(img_mask, img_beta1, img_beta2);
                varargout{1} = BOLD_beta;

            case 'activationVolume'
                activatedMap = spm_vol([mainStruct.meta.folder '\' nam '\derived\res\spmT_0001.nii']);
                makeSameResolution(activatedMap.fname, Vox_mask.fname);
                Vox_mask = spm_vol([mainStruct.meta.folder '\' nam '\anat\derived\rsw' nam '_sp_mask_sham.nii']);
                actMap = spm_read_vols(activatedMap); 
                %actMap(actMap< 3.145452) = 0; actMap(actMap>= 3.145452) = 1; % 3.145452  - from SPM results page, it is slightly different across subjexts
                Vox_map = spm_read_vols(Vox_mask);
                intersectionCoef = sum(Vox_map(actMap>= 3.145452), 'all')/sum(Vox_map, 'all');
                varargout{1} = intersectionCoef;
                mainStruct.(nam).proc.bold.ActivationVolume = intersectionCoef;
                mainStruct = hp_make('save', mainStruct);



        end

    %% Get data from the  DATABASE

    case 'getValue'
        %use as [~, resTable] = hp_make('getValue', id, valueChain)
        %e.g. hp_make('getValue', 20, {'proc','sham', 'all', 'Glx_AC'});
        % gives you Glx absolute concentration value for 20th subject
        mainStruct = hp_make('load');
        id = varargin{1};
        valueChain = varargin{2};
        nam = sprintf('sub_%02i', id);
        
        temp_struct = mainStruct.(nam);
        for i=1:length(valueChain)
            if ~isfield(temp_struct, valueChain{i})
%                 error("There is no such value");
                return;
            else
                temp_struct = temp_struct.(valueChain{i});
            end
        end
        varargout{1} = temp_struct;

        %concentration quantification
    case 'AbsoluteConc'
        %use as [~, CONC ] = hp_make('AbsoluteConc', id,  condition, met)
        %formula according DOI: 10.1007/s10334-012-0305-z
        %formula (2) p.333
        mainStruct = hp_make('load');
        id = varargin{1};
        condition = varargin{2};
        met = varargin{3};
        nam = sprintf('sub_%02i', id);

        LCmodelConc = mainStruct.(nam).proc.(condition).all.(met);
        met_pars{1} = met;
        met_pars{2} = mainStruct.(nam).proc.(condition).all.GM;
        met_pars{3} = mainStruct.(nam).proc.(condition).all.WM;
        met_pars{4} = mainStruct.(nam).proc.(condition).all.CSF;

        CONC = AbsoluteConcQunatification(LCmodelConc, met_pars);
        mainStruct.(nam).proc.(condition).all.([met '_AC']) = CONC;
        mainStruct.(nam).proc_check.absConc = 1;
        varargout{1}= CONC;

        mainStruct = hp_make('save', mainStruct);

%         IntRatio = I*N1*1*1/(2*35880*0.7);
%         CONC = IntRatio*(2/N_H)*CONC_wat*(frac_GM*R_GM*CONT_GM+frac_WM*R_WM*CONT_WM+frac_CSF*R_CSF*CONT_CSF)/...
%             (met_frac_GM*met_R_GM+frac_WM*R_WM);

    case 'importResultsData2'
        %use as [~, resTable] = hp_make('importResultsData2', ResCase)
        %This case uses new function to get values from meta-struture
        mainStruct = hp_make('load');
        ResCase = varargin{1};
        switch ResCase
            case 'spectra_ALL'
                %This case get data for ALL spectral data comparison (without dynamics)
                Values = {'NAA', 'Cr', 'Glx' ,'NAAerr', 'Crerr','Glxerr' ,'LWCr' ,'LWNAA','GM' ,'WM' ,'CSF'};
                for i=4:32
                    for ii=1:length(Values)
                        valueChain = {'proc','act', 'all', Values{ii} };
                        tableColumns{2*ii-1, 1} = ['act_' Values{ii}];
                        [~, resTable(i, 2*ii-1)] = hp_make('getValue', i, valueChain);
                        valueChain = {'proc','sham', 'all', Values{ii}};
                        tableColumns{2*ii, 1} = ['sham_' Values{ii}];
                        [~, resTable(i, 2*ii)] = hp_make('getValue', i, valueChain);
                    end
                end
                resTable = array2table(resTable, 'VariableNames', tableColumns);
                varargout{1} = resTable;

                 writetable(resTable, 'E:\Alex\fMRS-heatPain\_meta\spectraAll.csv');
            case 'spectra_dynamics'
                %This case get data for spectral data comparison in
                %dynamics
                Values = {'NAA', 'Cr', 'Glx','NAAerr', 'Crerr','Glxerr' };
                for i=4:32
                    k = 1;
                    for ii=1:length(Values)
                        for ij = 1:6
                            sp_nam = sprintf('tp_%02i', ij);
                            valueChain = {'proc','act', sp_nam, Values{ii} };
                            tableColumns{1, k} = ['act_' Values{ii} '_' sp_nam];
                            [~, resTable(i, k)] = hp_make('getValue', i, valueChain);
                            k = k+1;
                            valueChain = {'proc','sham', sp_nam, Values{ii}};
                            tableColumns{1, k} = ['sham_' Values{ii} '_' sp_nam];
                            [~, resTable(i, k)] = hp_make('getValue', i, valueChain);
                            k = k+1;
                        end
                    end
                end
                resTable = array2table(resTable, 'VariableNames', tableColumns);
                varargout{1} = resTable;

%                 writetable(resTable, 'C:\Users\Science\YandexDisk\Work\data\fMRS-hp\results\spectraDynamic_sm.csv');

            case 'spectra_6p'
                Values = {'NAA', 'Cr', 'Glx' ,'NAAerr', 'Crerr','Glxerr' ,'LWCr' ,'LWNAA'};
                for i=4:32
                    k=1;
                    for ii=1:length(Values)
                        
                        for ij = 1:3
                            sp_nam = sprintf('p%i', ij);
                        
                            valueChain = {'proc','act', sp_nam, Values{ii} };
                            tableColumns{k, 1} = ['act_' sp_nam '_' Values{ii}];
                            [~, resTable(i, k)] = hp_make('getValue', i, valueChain);
                            k=k+1;
                            valueChain = {'proc','sham', sp_nam, Values{ii}};
                            tableColumns{k, 1} = ['sham_' sp_nam '_' Values{ii}];
                            [~, resTable(i, k)] = hp_make('getValue', i, valueChain);
                            k=k+1;
                        end
                    end
                end
                resTable = array2table(resTable, 'VariableNames', tableColumns);
                varargout{1} = resTable;

                 writetable(resTable, 'C:\Users\Science\YandexDisk\Work\data\fMRS-hp\results\spectra6P.csv');

            case 'reactionMRS'
                Values = {'rating', 'reaction_time'};
                for i=4:32
                    k=1;
                    ii = 1;
                    valueChain = {'proc','MRI', Values{ii} };
                    tableColumns{k, 1} =  Values{ii};
                    try
                        [~, rating] = hp_make('getValue', i, valueChain);
                    catch ME
                        if (strcmp(ME.message,'There is no such value'))
                            continue
                        end
                    end 
                    
                    resTable(i, k) = mean(rating);
                    k=k+1;
                    ii = 2;
                    valueChain = {'proc','MRS', Values{ii} };
                    tableColumns{k, 1} =  'meanReactTime';
                    
                    try
                        [~, reactTime] = hp_make('getValue', i, valueChain);
                    catch ME
                        if (strcmp(ME.message,'There is no such value'))
                            continue
                        end
                    end

                    resTable(i, k) = mean(reactTime);
                    k=k+1;
                    tableColumns{k, 1} =  'STD_ReactTime';
                    resTable(i, k) = std(reactTime);
                    k=k+1;
                    tableColumns{k, 1} =  'N_ReactTime';
                    resTable(i, k) = length(reactTime);
                    k=k+1;

                end
                resTable = array2table(resTable, 'VariableNames', tableColumns);
                varargout{1} = resTable;
                writetable(resTable, 'E:\Alex\fMRS-heatPain\_meta\est_MRI.csv');
                
               
            case 'BOLD_MRS'
                % only linewidthes in both cases for time points
                % [~, resTable] = hp_make('importResultsData2', 'BOLD_MRS');
                Values = {'LWCr' ,'LWNAA', 'HCr', 'HNAA', 'SNR_NAA'};
                for i=4:32
                    k=1;
                    for ii=1:length(Values)
                        for ij = 1:6
                            %Choose case carefully here
                            sp_nam = sprintf('tp_%02i', ij); %intial data
%                             sp_nam = sprintf('tp_%02i_sm', ij);%temporally smoothed data
                        
                            valueChain = {'proc','act', sp_nam, Values{ii} };
                            tableColumns{k, 1} = ['act_' sp_nam '_' Values{ii}];
                            [~, resTable(i, k)] = hp_make('getValue', i, valueChain);
                            k=k+1;
                            valueChain = {'proc','sham', sp_nam, Values{ii}};
                            tableColumns{k, 1} = ['sham_' sp_nam '_' Values{ii}];
                            [~, resTable(i, k)] = hp_make('getValue', i, valueChain);
                            k=k+1;
                        end
                    end
                end

                resTable = array2table(resTable, 'VariableNames', tableColumns);
                varargout{1} = resTable;

                writetable(resTable, 'C:\Users\Science\YandexDisk\Work\data\fMRS-hp\results\BOLD_MRS.csv');
                
            case 'BOLD_HRF_mrs'
                valueChain = {'proc','hfr_mrs'};
                for i=4:32
                    k=1;
                    for ij = 1:6
                        sp_nam = sprintf('%02i', ij);
                        tableColumns{k, 1} = [sp_nam];
                        k = k+1;
                    end
                        [~, resTable(i, :)] = hp_make('getValue', i, valueChain);
                end
                resTable = array2table(resTable, 'VariableNames', tableColumns);
                varargout{1} = resTable;
            case 'fMRI_metrics'
                Values = {'MeanIL','MeanIR','MeanSMA','ActivationVolume'};
                for i=4:32
                    k=1;
                    for ii = 1:4
                        valueChain = {'proc','bold' Values{ii} };
                        tableColumns{k, 1} = Values{ii};
                        [~, resTable(i, k)] = hp_make('getValue', i, valueChain);
                        k = k+1;
                    end
                end
                resTable = array2table(resTable, 'VariableNames', tableColumns);
                varargout{1} = resTable;

                writetable(resTable, 'C:\Users\Science\YandexDisk\Work\data\fMRS-hp\results\fMRI_metrics.csv');
                
        end

                
    case 'BOLD_correction'
        %here using NAA nad Cr signal parameters BOLD change detected
        %(should be narrowing of the lines) and specifically corrected
        %(using FID-A line broadening)
        %use as [] = hp_make('BOLD_correction', id)
        mainStruct = hp_make('load');
        id = varargin{1};
        nam = sprintf('sub_%02i', id);

        %Find linewidths
        for j= 1:6
            tp_nam = sprintf('tp_%02i', j);
            LW_met(j, 1) = mainStruct.(nam).proc.act.(tp_nam).LWCr;
            LW_NAA(j, 1) = mainStruct.(nam).proc.act.(tp_nam).LWNAA;
        end
        %find difference between TP LW and mean LW
        LW_Cr_n = LW_met-mean(LW_met); LW_NAA_n = LW_NAA-mean(LW_NAA);
        %make correction according to the difference
        %         for i=1:mainStruct.meta.subNumbers
        nam = sprintf('sub_%02i', id);
        if mainStruct.(nam).data_check.sp>0
            txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'a');
            fprintf(txt_protocol, '------ \n bold correction Cr signal: case - %s \n', nam);
            for ii=1:6
                tp_nam = sprintf('tp_%02i', ii);
                sp_name = [mainStruct.meta.folder '\' nam '\sp\derived\' nam '_all_' tp_nam '_act.SDAT'];
                sp = io_loadspec_sdat(sp_name, 1);
                sp = op_filter(sp, -LW_Cr_n(ii, 1));
                new_sp_name = [mainStruct.meta.folder '\' nam '\sp\derived\' nam '_all_' tp_nam '_act_bc.SDAT'];
                mrs_writeSDAT(new_sp_name, sp.fids);
                copyfile([sp_name(1:end-4) 'SPAR'], [new_sp_name(1:end-4) 'SPAR']);
                
                %write LW-changes into the protocol of mrs processing
                fprintf(txt_protocol, '%s \n', datetime("today"));
                fprintf(txt_protocol, 'linewidth change for tp_%02i: %f \n', ii, -LW_Cr_n(ii, 1));
                
                

            end
            % check new linewidth
            fclose(txt_protocol);
            hp_make('LinewidthAssessment', id, 'act', 1, 1);
            mainStruct = hp_make('load');
            for j= 1:6
                tp_nam = sprintf('tp_%02i', j);
                LW_met_bc(j, 1) = mainStruct.(nam).proc.act.(tp_nam).LWCr_bc;
                LW_NAA_bc(j, 1) = mainStruct.(nam).proc.act.(tp_nam).LWNAA_bc;
            end
            plot(1:6, [LW_met, LW_met_bc]);
            mainStruct.(nam).proc_check.bold_correction = 1;
            mainStruct = hp_make('save', mainStruct);
            
        end
                   




                
        case 'getResBOLD'
        %use as [~, resTable] = hp_make('getResBOLD')
        %gives results of the functional MRI experiment as a matrix
        mainStruct = hp_make('load');
        resTable = zeros(1,1);

        for id=1:mainStruct.meta.subNumbers
            nam = sprintf('sub_%02i', id);
            

           if isfield(mainStruct.(nam).proc, 'bold') 
              resTable(id, 1) = mainStruct.(nam).proc.bold.mean_delta_beta1;
           end
           
        end
        varargout{1} = resTable;



%% Quality control procedures
        % Quality control for fMRI data (EVoronkova 12122023)

    case 'Q6'

        % to activate the case -    hp_make('Q6',8)
        
        id = varargin{1};
        mainStruct = hp_make('load');
        nam = sprintf('sub_%02i', id);

        fmri_file = dir([mainStruct.meta.folder '\' nam '\derived\r' nam '_func.nii']);
        if isempty(fmri_file)
            fmri_file = dir([mainStruct.meta.folder '\' nam '\derived\ra' nam '_func.nii']);
        end
        rp_file = dir([mainStruct.meta.folder '\' nam '\derived\rp_' nam '_func.txt']);
        if isempty(rp_file)
            rp_file = dir([mainStruct.meta.folder '\' nam '\derived\rp_a' nam '_func.txt']);
        end
        SPM_file = dir([mainStruct.meta.folder '\' nam '\derived\res\SPM.mat']);

        fmri_file = fullfile(fmri_file(1).folder, fmri_file(1).name);
        rp_file = fullfile(rp_file(1).folder, rp_file(1).name);
        SPM_file = fullfile(SPM_file(1).folder, SPM_file(1).name);
        %fmri_file = ('C:\Helen\Docs\DHiT\Alex\FMRS\sub_08\func\sub_08_func.nii');
        %rp_file = ('C:\Helen\Docs\DHiT\Alex\FMRS\sub_08\derived\rp_sub_08_func.txt');
        %SPM_file = ('C:\Helen\Docs\DHiT\Alex\FMRS\sub_08\derived\res\SPM.mat');

        Q6_time_series_task(fmri_file, rp_file, SPM_file);
        saveas(gcf,[mainStruct.meta.folder '\' nam '\meta\Q6.jpg'],'jpg');


    case 'b1/b8'

%         hp_make('b1/b8',8)

        id = varargin{1};
        mainStruct = hp_make('load');
        nam = sprintf('sub_%02i', id);

        b1_file = ([mainStruct.meta.folder '\' nam '\derived\res\beta_0001.nii']);
        b8_file = ([mainStruct.meta.folder '\' nam '\derived\res\beta_0008.nii']);
        InsulaL_mask_file = ([mainStruct.meta.folder '\_meta\atlas_map\atlases_Lena\rInsula_left_thr20.nii']);
        InsulaR_mask_file = ([mainStruct.meta.folder '\_meta\atlas_map\atlases_Lena\rInsula_right_thr20.nii']);
        SMA_mask_file = ([mainStruct.meta.folder '\_meta\atlas_map\atlases_Lena\rSMA_thr20.nii']);

        HeaderInfoB1 = spm_vol(b1_file);
        %disp(HeaderInfoB1);
        HeaderInfoB8 = spm_vol(b8_file);
        HeaderInfoIL = spm_vol(InsulaL_mask_file);
        HeaderInfoIR = spm_vol(InsulaR_mask_file);
        HeaderInfoSMA = spm_vol(SMA_mask_file);
        %disp(HeaderInfoSMA);
        
        if HeaderInfoB1.dim ~= HeaderInfoIL.dim
            error('The dimensions of main file and mask_files are NOT equal')
        end

        DataB1 = spm_read_vols(HeaderInfoB1);
        DataB8 = spm_read_vols(HeaderInfoB8);
        DataIL = spm_read_vols(HeaderInfoIL);
        DataIR = spm_read_vols(HeaderInfoIR);
        DataSMA = spm_read_vols(HeaderInfoSMA);

        %Data_B1_B8 = DataB1./DataB8;
        Data_B1_B8 = zeros(size(DataB1));
        S_IL = 0;
        N_IL = 0;
        S_IR = 0;
        N_IR = 0;
        S_SMA = 0;
        N_SMA = 0;
        
        for i = 1:HeaderInfoB1.dim(1)
            for j = 1:HeaderInfoB1.dim(2)
                for k = 1:HeaderInfoB1.dim(3)
                    if DataIL(i,j,k)>0 & ~(isnan(DataB1(i,j,k))) & ~(isnan(DataB8(i,j,k)))
                        Data_B1_B8(i,j,k) = DataB1(i,j,k)/DataB8(i,j,k);
                        N_IL = N_IL+1;
                        S_IL = S_IL+Data_B1_B8(i,j,k);
                    end
                    if DataIR(i,j,k)>0 & ~(isnan(DataB1(i,j,k))) & ~(isnan(DataB8(i,j,k)))
                        Data_B1_B8(i,j,k) = DataB1(i,j,k)/DataB8(i,j,k);
                        N_IR = N_IR+1;
                        S_IR = S_IR+Data_B1_B8(i,j,k);
                    end
                    if DataSMA(i,j,k)>0 & ~(isnan(DataB1(i,j,k))) & ~(isnan(DataB8(i,j,k)))
                        Data_B1_B8(i,j,k) = DataB1(i,j,k)/DataB8(i,j,k);
                        N_SMA = N_SMA+1;
                        S_SMA = S_SMA+Data_B1_B8(i,j,k);
                    end
                end
            end
        end
        MeanIL = S_IL/N_IL; mainStruct.(nam).proc.bold.MeanIL = MeanIL;
        MeanIR = S_IR/N_IR; mainStruct.(nam).proc.bold.MeanIR = MeanIR;
        MeanSMA = S_SMA/N_SMA; mainStruct.(nam).proc.bold.MeanSMA = MeanSMA;
        %fprintf('IL: %6.4f, IR: %6.4f, SMA:%6.4f,\n', MeanIL, MeanIR, MeanSMA);
        fprintf('%6.4f\t %6.4f\t %6.4f\n', MeanIL, MeanIR, MeanSMA);

%     case 'all_b1/b8'
%         for i = 3:25
%             hp_make('b1/b8',i)
%         end
        hp_make('save', mainStruct);

    case 'MaskMeanSignal'
        % Insula (R ans L) and supplementary motor area (SMA)

        id = varargin{1};
        mainStruct = hp_make('load');
        nam = sprintf('sub_%02i', id);


        wr_func_file = ([mainStruct.meta.folder '\' nam '\derived\wr' nam '_func.nii']);
        InsulaL_mask_file = ([mainStruct.meta.folder '\_meta\atlas_map\atlases_Lena\rInsula_left_thr20.nii']);
        InsulaR_mask_file = ([mainStruct.meta.folder '\_meta\atlas_map\atlases_Lena\rInsula_right_thr20.nii']);
        SMA_mask_file = ([mainStruct.meta.folder '\_meta\atlas_map\atlases_Lena\rSMA_thr20.nii']);

        wr_func_file = ([mainStruct.meta.folder '\FMRS\' nam '\derived\wr' nam '_func.nii']);
        InsulaL_mask_file = ([mainStruct.meta.folder '\FMRS\Atlases\rInsula_left_thr20.nii']);
        InsulaR_mask_file = ([mainStruct.meta.folder '\FMRS\Atlases\rInsula_right_thr20.nii']);
        SMA_mask_file = ([mainStruct.meta.folder '\FMRS\Atlases\rSMA_thr20.nii']);


        HeaderInfo0 = spm_vol(wr_func_file);
        disp(HeaderInfo0(1));
        HeaderInfoIL = spm_vol(InsulaL_mask_file);
        %disp(HeaderInfoIL);
        HeaderInfoIR = spm_vol(InsulaR_mask_file);
        %disp(HeaderInfoIR);
        HeaderInfoSMA = spm_vol(SMA_mask_file);
        %disp(HeaderInfoSMA);

        if HeaderInfo0(1).dim ~= HeaderInfoIL.dim
            error('The dimensions of the wr_func_file and mask_files are NOT equal')
        end

        Data0 = spm_read_vols(HeaderInfo0);
        DataIL = spm_read_vols(HeaderInfoIL);
        DataIR = spm_read_vols(HeaderInfoIR);
        DataSMA = spm_read_vols(HeaderInfoSMA);

        L = length(Data0);
        S_IL = zeros(1,L);
        N_IL = 0;
        S_IR = zeros(1,L);
        N_IR = 0;
        S_SMA = zeros(1,L);
        N_SMA = 0;
        for i = 1:HeaderInfo0(1).dim(1)
            for j = 1:HeaderInfo0(1).dim(2)
                for k = 1:HeaderInfo0(1).dim(3)
                    if DataIL(i,j,k)>0
                        N_IL = N_IL+1;
                        for x = 1:L
                            S_IL(x) = S_IL(x)+Data0(i,j,k,x);
                        end
                    end
                    if DataIR(i,j,k)>0
                        N_IR = N_IR+1;
                        for x = 1:L
                            S_IR(x) = S_IR(x)+Data0(i,j,k,x);
                        end
                    end
                    if DataSMA(i,j,k)>0
                        N_SMA = N_SMA+1;
                        for x = 1:L
                            S_SMA(x) = S_SMA(x)+Data0(i,j,k,x);
                        end
                    end
                end
            end
        end
        MeanIL = S_IL/N_IL;
        writematrix(MeanIL,[mainStruct.meta.folder '\FMRS\' nam '\derived\MeanIL.txt'],"Delimiter",'tab');
        MeanIR = S_IR/N_IR;
        writematrix(MeanIR,[mainStruct.meta.folder '\FMRS\' nam '\derived\MeanIR.txt'],"Delimiter",'tab');
        MeanSMA = S_SMA/N_SMA;
        writematrix(MeanSMA,[mainStruct.meta.folder '\FMRS\' nam '\derived\MeanSMA.txt'],"Delimiter",'tab');

    case 'HRF_GLMsingle'
        % to activate the case -    hp_make('HRF_GLMsingle',12)
        % git clone --recurse-submodules https://github.com/cvnlab/GLMsingle.git          
        id = varargin{1};
        mainStruct = hp_make('load');
        nam = sprintf('sub_%02i', id);

        timingInfo = heatPain_makeRegressor([mainStruct.meta.folder '\' nam '\func\' nam '_bold.xlsx']);
        startTime = getTTLtime([mainStruct.meta.folder '\' nam '\func\' nam '_bold.xlsx']);
        timingInfo(:, 1) = timingInfo(:, 1) - startTime - (12 - mainStruct.(nam).proc.dummy_time);

        %timingInfo(:,1) %stimulus start time
        %timingInfo(:,2) %stimulus duration
        stimdur = mean(timingInfo(:,2));

        wr_func_file = ([mainStruct.meta.folder '\' nam '\derived\r' nam '_func.nii']);
        HeaderInfo0 = spm_vol(wr_func_file);
        data0 = spm_read_vols(HeaderInfo0);
        tr = HeaderInfo0(1).private.timing.tspace;
        
        design1 = round(timingInfo(:,1));
        design_tr1 = fix(design1/tr)+1;
        design_tr = design_tr1(design_tr1 <= size(data0, 4));
        S = sparse(size(data0, 4),2);
        for i = design_tr
            S(i,1) = 1;   
        end

%         for i = design_tr
%             if i <= size(data, 4)
%                 S(i,1) = 1;
%             end
%         end

        design = cell(1);
        design{1,1} = S;
        %disp(design{1,1});
        data = cell(1);
        data{1,1} = data0;

        if ~exist([mainStruct.meta.folder '\' nam '\GLMsingle'])
            mkdir ([mainStruct.meta.folder '\' nam '\GLMsingle'])
        end
        outputdir = ([mainStruct.meta.folder '\' nam '\GLMsingle']);

        %opt = struct('wantmemoryoutputs',[0 1 0 1]); % by default 0001
        opt.wantglmdenoise = 0;     % by default 1
        opt.wantfracridge = 0;      % by default 1
        opt.wantmemoryoutputs = [0 1 0 0]; % by default 0001
        opt.wantfileoutputs = [1 1 0 0]; %  by default 1111
     

        [results] = GLMestimatesingletrial(design,data,stimdur,tr,outputdir,opt);
        
    case '2lvlanalysis'
   % hp_make('2lvlanalysis', temp, hand);
   mainStruct = hp_make('load');
   if length(varargin)>0
       temp = varargin{1};
       if length(varargin)<3
           hand = varargin{2};
       else
           other = varargin{3};
       end
   end
    %special case (29-07) for estimation of the stimuli
    
    if true
        %get from csv file
        estTab = readtable([mainStruct.meta.folder '\_meta\est_MRI.csv']);
        estTab = table2array(estTab);
        estTab(isnan(estTab)) = 0;
        regrTab = []; k =1;
        for i=3:32
            %find those data that dont have null
            nam = sprintf('sub_%02i', i);
            estTab(isnan(estTab)) = 0;
            if ~isempty(find(estTab(i, :), 1))
                input{1, 1}{k, 1} = [mainStruct.meta.folder '\' nam '\derived\res\con_0001.nii']; 
                regrTab = [regrTab; estTab(i, :)];
                Temp(k, 1) = mainStruct.(nam).proc.selected_temp;k = k+1;
            end
        end
        input{1, 2} = regrTab(:, 2);
    end
   for i=3:32
       nam = sprintf('sub_%02i', i);
       input{1, 1}{i-2, 1} = [mainStruct.meta.folder '\' nam '\derived\res\con_0001.nii'];
   end
   input{1, 2} = temp;
   input{1, 3} = hand;

   secondLvl(input);
   




        %% Some not very important features or minor experiments (consider to remove it)
    case 'GLM_estimation'
        %Here in GLM design start of the stimulus presentation is engaged
        %with button response
        hp_default;
        id = varargin{1};
        nam = sprintf('sub_%02i', id);
        mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\derived\res_exp']);
        nrun = 1; % enter the number of runs here
        jobfile = {[mainStruct.meta.folder '\_meta\stats_empty_job.m']};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(5, nrun);
        for crun = 1:nrun
            inputs{1, crun} = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\res_exp']}; % fMRI model specification: Directory - cfg_files
            fmri_img = spm_vol([mainStruct.meta.folder mainStruct.(nam).folder '\derived\swra' nam '_func.nii']); NSA = length(fmri_img);
            for i=1:NSA
                func_data(i, 1) = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\swra' nam '_func.nii,' num2str(i)]};
            end
            inputs{2, crun} = func_data; % fMRI model specification: Scans - cfg_files
            if floor(mod(mainStruct.(nam).data_check.funcTable, 100)/10)>0
                startTime = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
                regressorList = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
                if length(mainStruct.(nam).proc.react_time)<length(regressorList)
                    delta = length(regressorList)-length(mainStruct.(nam).proc.react_time);
                    mainStruct.(nam).proc.react_time(length(mainStruct.(nam).proc.react_time)+1:length(mainStruct.(nam).proc.react_time)+delta) = mean(mainStruct.(nam).proc.react_time);
                    regressorList(:, 1) = startTime - startTime(1, 1) - (12 - mainStruct.(nam).proc.dummy_time)+mainStruct.(nam).proc.react_time;
                else
                    regressorList(:, 1) = startTime - startTime(1, 1) - (12 - mainStruct.(nam).proc.dummy_time)+mainStruct.(nam).proc.react_time;
                end
            end
            inputs{3, crun} = regressorList(:, 1); % fMRI model specification: Onsets - cfg_entry
            inputs{4, crun} = regressorList(:, 2); % fMRI model specification: Durations - cfg_entry
            inputs{5, crun} = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\rp_a' nam '_func.txt']}; % fMRI model specification: Multiple regressors - cfg_files
        end
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs, inputs{:});

    case 'BOLD effect estimate'
        %Here BOLD-effect using MRS daa was estimated.
        % these five lines is dedicated to obtain linewidth and normilise
        % it as: (X-mean(X))/std(X); (mean and std by rows)
        cases = {'sham', 'act'};
        mets = {'NAA', 'Cr'};
        for i=1:2
            [~, resTable] = hp_make('getResSP', cases{i}, 'LW');
            for ii=1:2
                LW_met = resTable.LW.(mets{ii});
                LW_met = LW_met(5:18, :);
%                 LW_met(7, :) = [];
                LW_met_norm = (LW_met-repmat(mean(LW_met, 2), 1, 6))./repmat(std(LW_met, [], 2), 1, 6);
                LW_case_met{i, ii} = LW_met_norm;
            end
        end
        
        %This is theroretically divided BOLD (got by standard-hrf convolved with
        %stimulus funcion)
%         subjs = [5:10, 12:18];
        subjs = [5:18];
        for i=1:13
            nam = sprintf('sub_%02i', subjs(i));
            hrf_mrs(:, i) = mainStruct.(nam).proc.hfr_mrs;
        end
        hrf_mrs = hrf_mrs';

        varargout{1} = LW_case_met;


    case 'expDynamicsSNR'
        mainStruct = hp_make('load');
        pathToSave = 'G:\_other\fMRS-heatPain\_unsorted\dyn60\';
        %open only SHAM data, preprocess it and save by 60 averaged dynamics
        for i=3:mainStruct.meta.subNumbers
            nam = sprintf('sub_%02i', i);
            sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\' nam '_sham.SDAT'], 1);
            pars{1, 1} = [1.9 2.1 1]; 
            sp_fpa_1 = op_freqAlignAverages_fd(sp, pars{1, 1}(1), pars{1, 1}(2), pars{1, 1}(3), 'n');
            for ii = 3:6
                sp_fpa_dyn_1{ii,1} = op_takeaverages(sp, [ii:4:240]);
                sp_fpa_av{ii,1} = op_averaging(sp_fpa_dyn_1{ii,1});
                pars{2, 1} = {[4.4 5], 20};
                sp_fpa_av_wr{ii,1} = op_removeWater(sp_fpa_av{ii,1}, pars{2, 1}{1}, pars{2, 1}{2});
                mrs_writeSDAT([pathToSave nam '_' num2str(ii) '.SDAT'],...
                    sp_fpa_av_wr{ii,1}.fids);
                copyfile([mainStruct.meta.folder '\' nam '\sp\' nam '_sham.SPAR'],...
                    [pathToSave nam '_' num2str(ii) '.SPAR']);
                my_editSPAR([pathToSave nam '_' num2str(ii) '.SPAR'],{'rows','spec_num_row', 'spec_row_upper_val'}, {1,1,1});
            end
        end

    case 'sub11_waterBOLD'
        % SPECIAL_CASE: subject 11 has unsupressed spectra made with heat
        % pain stimulation (only act state, sham is OK)
        %here linewidth and height of the water signal (4.65) was found
        % use as [LW H] = hp_make('sub11_waterBOLD');

        mainStruct = hp_make('load');
        nam = sprintf('sub_%02i', 11);
        sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\' nam '_act.SDAT'], 1);
%         for i=1:6
%             sp_nam = sprintf('tp_%02i', i);
%             sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\derived\' nam '_all_' sp_nam '_act.SDAT'], 1);
%             LW(i,1) = op_getLW(sp, 4.2, 5, 16);
%             HGHT(i, 1) = op_getPeakHeight(sp, 4.2, 5);
%         end
        for i=1:sp.sz(2)
            sp1= op_takeaverages(sp, i);
            LW(i,1) = op_getLW(sp1, 4.2, 5, 16);
            HGHT(i, 1) = op_getPeakHeight(sp1, 4.2, 5);
        end
        varargout{1} = LW;
        varargout{2} = HGHT;

        regressorList = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_fmrs.xlsx']);
        regressorList(:, 1) = regressorList(:, 1) - 17.87;
        load('F:\_other\fMRS-heatPain\sub_11\spec_case\SPM.mat');
        HRF = SPM.xX.X(:,1); 

        resp_norm = LW/mean(LW)-1; resp_norm = circshift(resp_norm, 1);
        resp_norm2 = H/mean(H)-1;resp_norm2 = circshift(resp_norm2, 1);

        [r, p] = corr(HRF, resp_norm);
        [r2, p2] = corr(HRF, resp_norm2);

    case 'VASassessment'
        % [~, VAStable] = hp_make('VASassessment')
        %Here VAS estimates are assessed from pretest files and using
        %tempereature_approximation VAS estimate was obtained
        mainStruct = hp_make('load');

        for i=1:mainStruct.meta.subNumbers
            nam = sprintf('sub_%02i', i);
            %already done
%             hp_make('new_field', 'mainStruct.(nam).proc.estimate', 7);
%             hp_make('new_field', 'mainStruct.(nam).proc_check.estimates', 0);
%             hp_make('new_field', 'mainStruct.(nam).proc.VAS', 7);
            switch i
                case 5
                    mainStruct.(nam).proc.temp = [39 42 45 42 39 43 37 44 42];
                    mainStruct.(nam).proc.estimate = [5 6 8.5 5 3 7 2 8 8];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 43;
                case 6
                    mainStruct.(nam).proc.temp = [39 42 45 46 47 48];
                    mainStruct.(nam).proc.estimate = [2 5 6 6 7 8.5];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 47;
                case 7
                    mainStruct.(nam).proc.temp = [39 42 45 46 47 48 42 48 50 49];
                    mainStruct.(nam).proc.estimate = [1 2 3 3 4 6 1 6 8 8];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 48;
                case 8
                    mainStruct.(nam).proc.temp = [39 42 45 46 47 48 41 39 47 46 49 50];
                    mainStruct.(nam).proc.estimate = [1, 3, 3 4 5 6 1 0 4 4 7 8];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 49;
                case 10
                    mainStruct.(nam).proc.temp = [39 42 45 46 47];
                    mainStruct.(nam).proc.estimate = [2 4 7 8 9];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 46;
                case 12
                    mainStruct.(nam).proc.temp = [39 42 45 46 47];
                    mainStruct.(nam).proc.estimate = [5 6 7 7 8];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 48;
                case 14
                    mainStruct.(nam).proc.temp = [39 42 45 46 47 48 41 45 49 48];
                    mainStruct.(nam).proc.estimate = [3 4 4 5 5 6 3 4 6 7.5];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 48;
                case 15
                    mainStruct.(nam).proc.temp = [37 39 41 43 45 47];
                    mainStruct.(nam).proc.estimate = [0 1 3 4 5 7];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 47;
                case 16
                    mainStruct.(nam).proc.temp = [39 42 45 46 47 48 41 43 47 45 49 50];
                    mainStruct.(nam).proc.estimate = [1 2 3 4 4 5 2 3 4 4 5 6];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 50;
                case 17
                    mainStruct.(nam).proc.temp = [39 42 45 46 47];
                    mainStruct.(nam).proc.estimate = [3 5 6 7 8];
                    mainStruct.(nam).proc_check.estimates = 1;
                    mainStruct.(nam).proc.selected_temp = 47;

                end
        end
        for i=1:mainStruct.meta.subNumbers
            nam = sprintf('sub_%02i', i);
            if mainStruct.(nam).proc_check.estimates
                mainStruct.(nam).proc.VAS = temperature_approximation(mainStruct.(nam).proc.temp, mainStruct.(nam).proc.estimate, mainStruct.(nam).proc.selected_temp);
                VAStable(i).name = nam;
                VAStable(i).value = mainStruct.(nam).proc.VAS;
            end
        end
        varargout{1} = VAStable;
        hp_make('save', mainStruct);

    case 'BOLDextraction'
        %Here BOLD from three regions exctracted data was resampled
        mainStruct = hp_make('load');
        dats = readtable('C:\Users\Science\YandexDisk\Work\data\fMRS-hp\BOLD_mask.xlsx');
        dats(2,:) = [];
        for i=3:13
            bold_dat{i-2}=dats(i*3-7:i*3-5, 3:172);
            bold_dat{i-2} = table2array(bold_dat{i-2});
            bold_dat{i-2}(isnan(bold_dat{i-2})) = [];
            if size(bold_dat{i-2}, 1)<2
                bold_dat{i-2}= reshape(bold_dat{i-2}, 3, length(bold_dat{i-2})/3);
            end
        end
        for i=1:13-2
            x_len = length(bold_dat{i});
            for ii = 1:3
                bold_dat_rsmp{i}(ii, :) = interp1 ([0:3:(x_len-1)*3], bold_dat{i}(ii, :), [0:2:(x_len-1)*3], 'pchip');
                
            end
            bold_dat_rsmp{i}(4, :) = [0:2:(x_len-1)*3];
        end
        
        for i=1:11
            nam = sprintf('sub_%02i', i+2);
            av_resp(i, :) = zeros(1, 10);
            timingInfo = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
            startTime = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
%             switch i
%                 case 3
%                     startTime = 14.50;  %for some reason there is no tag in xlxs file
%                 case 4
%                     startTime = 28.024;
%                 otherwise
%                     startTime = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
%             end
            timingInfo(:, 1) = timingInfo(:, 1) - startTime - (12 - mainStruct.(nam).proc.dummy_time);
            timingInfo(:, 1) = round(timingInfo(:,1)/2);
            time_slice = timingInfo(timingInfo(:,1)<max(bold_dat_rsmp{i}(4, :))/2, 1);
            dyns_perstim = diff(time_slice);dyns_perstim(length(time_slice)) = max(bold_dat_rsmp{i}(4, :))/2-time_slice(end);
            x = zeros(10, 1);
            for ii=1:length(time_slice)                
                for ij=1:dyns_perstim(ii)
                    if  time_slice(ii)+ij-1>0
                        if ij>10
                            break
                        end
                        av_resp(i, ij) = av_resp(i, ij) + bold_dat_rsmp{i}(3, time_slice(ii)+ij-1);
                        x(ij, 1) = x(ij, 1)+1;
                    end
                end
                
            end
            for ij = 1:10
                if x(ij)<15
                    x(ij) = 0;
                    av_resp(i, ij) = 0;
                end
            end
            av_resp(i, x~=0) =av_resp(i, x~=0)./x(x~=0)';
        end
        for i=1:size(av_resp, 1)
            av_resp(i, :) = av_resp(i, :)./mean(av_resp(i, av_resp(i, :)~=0));
        end
        varargout{1} = av_resp;
        %after that
        % BOLD_response{1,1} = av_table;

    case 'meanStimDuaration'
        mainStruct = hp_make('load');
        for i=1:11
            nam = sprintf('sub_%02i', i+2);
            timingInfo = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
            startTime = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
            timingInfo(:, 1) = timingInfo(:, 1) - startTime - (12 - mainStruct.(nam).proc.dummy_time);
            stim_duartion(i,1) = mean(timingInfo(:,2));
        end
        varargout{1} = stim_duartion;

    case 'quantifyDiceCoeff'
        %[~, Dice] = hp_make('quantifyDiceCoeff', id)
        %Dice coefficient quantification between spectra localization in
        %ACT and REST conditions
        id = varargin{1};
        condition = {'act', 'sham'};
        mainStruct = hp_make('load');
        nam = sprintf('sub_%02i', id);
        
        if (mainStruct.(nam).proc_check.all.(condition{1})+mainStruct.(nam).proc_check.all.(condition{2}))>1
            img1 = spm_vol(mainStruct.(nam).proc.(condition{1}).all.path);
            img1 = spm_read_vols(img1);
            img2 = spm_vol(mainStruct.(nam).proc.(condition{2}).all.path);
            img2 = spm_read_vols(img2);
            mainStruct.(nam).proc.Dice = sum(2*img1.*img2, 'all')/sum(img1+img2, 'all');
        end
        varargout{1}=  mainStruct.(nam).proc.Dice;
        mainStruct = hp_make('save', mainStruct);

    case 'Test: pain responce'
        mainStruct = hp_make('load');

        
        estTab = readtable([mainStruct.meta.folder '\_meta\est_MRI.csv']);
        estTab = readtable([mainStruct.meta.folder '\_meta\est_MRI.csv']);
        estTab = table2array(estTab);
        estTab(isnan(estTab)) = 0;
        regrTab = []; k =1;
        for i=14:32
            %find those data that dont have null
            nam = sprintf('sub_%02i', i);
            fmri_img = spm_vol([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_func.nii']);
            NSA = length(fmri_img);
            if ~isempty(find(estTab(i, :), 1))
                mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\derived\res_test\']);
                nrun = 1; % enter the number of runs here
                jobfile = {[mainStruct.meta.folder '\_meta\stats_empty_job.m']};
                jobs = repmat(jobfile, 1, nrun);
                inputs = cell(5, nrun);
                for crun = 1:nrun
                    inputs{1, crun} = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\res_test']}; % fMRI model specification: Directory - cfg_files
                    for ij=1:NSA
                        func_data(ij, 1) = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\swra' nam '_func.nii,' num2str(i)]};
                    end
                    inputs{2, crun} = func_data; % fMRI model specification: Scans - cfg_files
                    if floor(mod(mainStruct.(nam).data_check.funcTable, 100)/10)>0
                        regressorList = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
                        %                 regressorList(:, 1) = regressorList(:, 1) - regressorList(1, 1) - (12 - mainStruct.(nam).proc.dummy_time);
                        %TTL_time = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
                        regressorList(:, 1) = mainStruct.(nam).proc.MRI.TTLtimes(1:length(regressorList))-mainStruct.(nam).proc.MRI.TTLtimes(1)+estTab(i, 2);
                    end
                    inputs{3, crun} = regressorList(:, 1); % fMRI model specification: Onsets - cfg_entry
                    inputs{4, crun} = regressorList(:, 2); % fMRI model specification: Durations - cfg_entry
                    inputs{5, crun} = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\rp_a' nam '_func.txt']}; % fMRI model specification: Multiple regressors - cfg_files
                end
                spm('defaults', 'FMRI');
                spm_jobman('run', jobs, inputs{:});

            end
        end

        
        
        



    case 'BugFix_LWbc_toBC'
        %Made 18-03-2024 AYakovlev
        %transfer LC values to different place in the structure
        mainStruct = hp_make('load');
        
        for i=1:mainStruct.meta.subNumbers
            nam = sprintf('sub_%02i', i);
            if mainStruct.(nam).proc_check.bold_correction
            for j=1:6
                tp_nam = sprintf('tp_%02i', j);
                mainStruct.(nam).proc.act.([tp_nam '_bc']).LWCr = mainStruct.(nam).proc.act.(tp_nam).LWCr_bc;
                mainStruct.(nam).proc.act.([tp_nam '_bc']).LWNAA = mainStruct.(nam).proc.act.(tp_nam).LWNAA_bc;
                mainStruct.(nam).proc.act.(tp_nam).LWCr_bc = [];
                mainStruct.(nam).proc.act.(tp_nam).LWNAA_bc = [];
            end
            end
        end
        hp_make('save', mainStruct);

         case 'EstimateIncreaseTime'
        %hp_make('EstimateIncreaseTime')
        mainStruct = hp_make('load');
        incTimeTemp = zeros(32, 2);
        for id=4:32
            
            nam = sprintf('sub_%02i', id);

            TSAtab = readtable([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
            temp= TSAtab.Tec_C_;
            time = TSAtab.Timestamp_msec_;

            time = time/1000;
            baselineTemp = median(temp(1:1000));
            maxTemp = quantile(temp,0.9);
            startInc = find(temp-baselineTemp>baselineTemp*0.05);startInc = time(startInc(1));
            endInc = find(maxTemp-temp<0);endInc = time(endInc(1));
            incTimeTemp(id, 1) = endInc - startInc; incTimeTemp(id, 2) = maxTemp;
        end
        varargout{1} = incTimeTemp;
    case 'setTemperature'
        %hp_make('setTemperature')
        mainStruct = hp_make('load');
        tab = readtable('E:\Alex\fMRS-heatPain\_meta\ResponsesAnalisys.xlsx');
        for i=1:32
            nam = sprintf('sub_%02i', i);
            if isnan(tab.Temp(i))
                continue;
            end
            mainStruct.(nam).proc.selected_temp = tab.Temp(i);
        end
         hp_make('save', mainStruct);

end
end

function Estimates = MRSassessment(sp)
sp.flags.averaged = 1;
sp.dims.averages = 0;
sp = op_autophase(sp, 1.8, 3.5, 0);

Estimates.LWCr = op_getLW(sp, 2.8, 3.1);
Estimates.LWNAA = op_getLW(sp, 1.8, 2.15);

Estimates.HCr = op_getPeakHeight(sp, 2.8, 3.1);
Estimates.HNAA = op_getPeakHeight(sp, 1.8, 2.15);


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
        thrs_1 = find(tstamp_1>0.2); thrs_1 = thrs_1(1);
        startTime = eventsTable.Timestamp_msec_(thrs_1)/1000;
    end
    

end

function callNormilise(mask_image, deformations_map)
    matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = {deformations_map};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample(1) = {mask_image};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

    spm_jobman('run',matlabbatch);
end

function callSmothing(img)
    matlabbatch{1}.spm.spatial.smooth.data(1) = {img};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

    spm_jobman('run',matlabbatch);
end

function makeSameResolution(img1, img2)
    matlabbatch{1}.spm.spatial.coreg.write.ref = {img1};
    matlabbatch{1}.spm.spatial.coreg.write.source = {img2};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
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


function callfMRIProcessing(inputs, slice_case, segment)
    mainStruct = hp_make('load');
    
    k=1;
    matlabbatch{k}.spm.temporal.st.scans = {inputs{1, 1}};
    matlabbatch{k}.spm.temporal.st.nslices = 35;
    matlabbatch{k}.spm.temporal.st.tr = 3;
    matlabbatch{k}.spm.temporal.st.ta = 2.91428571428571;
    switch slice_case
        case 'interleaved'
            matlabbatch{1}.spm.temporal.st.so = [1	7	13	19	25	31	2	8	14	20	26	32	3	9	15	21	27	33	4	10	16	22	28	34	5	11	17	23	29	35	6	12	18	24	30];
        case 'ascending'
            matlabbatch{1}.spm.temporal.st.so = [1:35];
    end
    matlabbatch{k}.spm.temporal.st.refslice = 1;
    matlabbatch{k}.spm.temporal.st.prefix = 'a';
    k=k+1;
    matlabbatch{k}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
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
    matlabbatch{k}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
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

function CONC = AbsoluteConcQunatification(LCmodelConc, met_pars)
    met = met_pars{1};
    frac_GM = met_pars{2};
    frac_WM = met_pars{3};
    frac_CSF = met_pars{4};

    switch met
        case 'Glx'
            N1 = 1;
            %accroding to https://www.mr.ethz.ch/abstracts/files/ismrm15_0202.pdf
            met_R_GM = exp(-35/144); %echo time/T2 of glu in GM
            met_R_WM = exp(-35/106); %echo time/T2 of glu in WM
            N_H = 1;
            met_frac_GM = 1;
            met_frac_WM = 1;
        case 'NAA'
            N1 = 1;
            %accroding to https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/mrm.21715
            met_R_GM = exp(-35/269); %echo time/T2 of NAA in GM
            met_R_WM = exp(-35/374); %echo time/T2 of NAA in WM
            N_H = 1;
            met_frac_GM = 1;
            met_frac_WM = 1;

        case 'Cr'
            N1 = 1;
            %accroding to https://www.mr.ethz.ch/abstracts/files/ismrm15_0202.pdf
            met_R_GM = exp(-35/156); %echo time/T2 of Cr in GM
            met_R_WM = exp(-35/179); %echo time/T2 of Cr in WM
            N_H = 1;
            met_frac_GM = 1;
            met_frac_WM = 1;
    end

    %Water parameters
    %https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/epdf/10.1002/nbm.4257
    %https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/mrm.20901
    CONC_wat = 53800;
    R2_GM = exp(-35/93); R1_GM = (1 - exp(-2000/832));  CONT_GM = 43300/CONC_wat;
    R2_WM = exp(-35/73); R1_WM = (1 - exp(-2000/1331)); CONT_WM = 36100/CONC_wat;
    R2_CSF = exp(-35/23); R1_CSF = (1 - exp(-2000/3817)); CONT_CSF = 53800/CONC_wat;

    


    IntRatio = LCmodelConc*N1*1*1/(2*35880*0.7); %denominator was taken from LCmodel manual
    CONC_csf_corrected = IntRatio/(1-frac_CSF)*CONC_wat;
    CONC = CONC_csf_corrected*(2/N_H)*(CONT_GM*frac_GM*R1_GM*R2_GM+frac_WM*R1_WM*R2_WM*CONT_WM+frac_CSF*R1_WM*R2_CSF*CONT_CSF)/...
        (met_frac_GM*met_R_GM+frac_WM*met_R_WM);
end

function [X] = makeHRF(onsets)

    %make correct onset times
%     [U] = spm_get_ons(SPM,s)
    pars = [6, 16, 1, 1, 6, 1, 32];
    %make basis function 
    global SPM1
    [bf, p]      = spm_hrf(2, pars(1:6), 16);
    bf = bf*pars(7);
    %convolve it with stimuli function from SPM struct (pre-defined using
    %GUI)  
    [X,Xn,Fc] = spm_Volterra(SPM1.Sess(1).U, SPM1.xBF.bf, SPM1.xBF.Volterra);
    %-Resample regressors at acquisition times (32 bin offset)
    %----------------------------------------------------------------------
    if ~isempty(X)
        X = X((0:(145 - 1))*SPM1.xBF.T + SPM1.xBF.T0 + 32,:);
    end
    X = X(:, 1);
end

function secondLvl(inputs)
matlabbatch{1}.spm.stats.factorial_design.dir = {'E:\Alex\fMRS-heatPain\_meta\2lvl\test29-07'};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = inputs{1, 1};
if size(inputs, 2)>1
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = inputs{1, 2};
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'temperature';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
end
if size(inputs, 2)>2
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = inputs{1, 3};
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'hands';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
end
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '1';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
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
matlabbatch{4}.spm.stats.results.export{2}.pdf = true;

spm_jobman('run',matlabbatch);
end