%This script is devoted to process all the data obtained with thermal pain
%stimuli

% AYakovlev 13-10-2023
% AYakovlev 13-11-2023

% id= 3;
% mainStruct = hp_make('sp_proc', id, 'sham', 1);
% mainStruct = hp_make('sp_proc', id, 'act', 1);
% hp_make('copyData', id,  'F:\heatPainSP', 'tp_mrs');
% hp_make('copyData', id,  'F:\heatPainSP', 'water_mrs');
% %%
% mainStruct = hp_make('load'); 
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
            eval([fieldname ' = ' num2str(default_value) ';']);
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
        if length(varargin)>0
            mainStruct = varargin{1};
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

        load([mainStruct.meta.folder '\_meta\spatial_proc.mat']); %loaded into matlabbatch
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
            jobfile = {'G:\_other\fMRS-heatPain\_meta\spatial_empty_job.m'};
            jobs = repmat(jobfile, 1, nrun);
            inputs = cell(3, nrun);
            for crun = 1:nrun
                inputs{1, crun} = func_data; % Realign: Estimate & Reslice: Session - cfg_files
                inputs{2, crun} = anat_data; % Coregister: Estimate: Reference Image - cfg_files
                inputs{3, crun} = anat_data; % Segment: Volumes - cfg_files
            end
            spm('defaults', 'FMRI');
            spm_jobman('run', jobs, inputs{:});


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
        jobfile = {'G:\_other\fMRS-heatPain\_meta\stats_empty_job.m'};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(5, nrun);
        for crun = 1:nrun
            inputs{1, crun} = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\res']}; % fMRI model specification: Directory - cfg_files
            for i=1:NSA
                func_data(i, 1) = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\swr' nam '_func.nii,' num2str(i)]};
            end
            inputs{2, crun} = func_data; % fMRI model specification: Scans - cfg_files
            if floor(mod(mainStruct.sub_03.data_check.funcTable, 100)/10)>0
                regressorList = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_bold.xlsx']);
                regressorList(:, 1) = regressorList(:, 1) - regressorList(1, 1) - (12 - mainStruct.sub_06.proc.dummy_time);
            end
            inputs{3, crun} = regressorList(:, 1); % fMRI model specification: Onsets - cfg_entry
            inputs{4, crun} = regressorList(:, 2); % fMRI model specification: Durations - cfg_entry
            inputs{5, crun} = {[mainStruct.meta.folder mainStruct.(nam).folder '\derived\rp_' nam '_func.txt']}; % fMRI model specification: Multiple regressors - cfg_files
        end
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs, inputs{:});

% MRS processing
        
    case 'sp_proc'
        %use as hp_make('sp_proc', id, sp_name, time_points_action)
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
        time_point_action = varargin{3};
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
        
        if time_point_action<1
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
            
        else
            %process as time points
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
                sp_fpa_av_wr = sp_out;

                % save intothe protocol
                txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'a');
                new_dir = [mainStruct.meta.folder '\' nam '\sp\derived\'];
                mkdir(new_dir);
                sp_naming = sprintf('%s_all_%s', nam, [tp_nam '_' sp_name]);
                mrs_writeSDAT([new_dir sp_naming '.SDAT'],...
                    sp_fpa_av_wr.fids);
                copyfile([mainStruct.meta.folder '\' nam '\sp\' nam '_' sp_name '.SPAR'],...
                    [new_dir sp_naming '.SPAR']);
                my_editSPAR([new_dir sp_naming '.SPAR'],{'rows','spec_num_row', 'spec_row_upper_val'}, {1,1,1});
                fprintf( txt_protocol, 'saved in: %s \n', [new_dir sp_naming '.SDAT']);
                fclose(txt_protocol);
            end
            mainStruct.(nam).proc_check.timepoints_spectra.(sp_name) = 1;
            mainStruct = hp_make('save', mainStruct);


           
        end
    
    case 'LinewidthAssessment'
        %use as hp_make('LinewidthAssessment', id, condition, tp)
        % assess FWHM of Cr and NAA signals
        mainStruct = hp_make('load');
        id = varargin{1};
        condition = varargin{2};
        tp = varargin{3};
        nam = sprintf('sub_%02i', id);
        
        if tp
            for i=1:mainStruct.(nam).proc_check.tp_spectra_res.(condition)
                sp_nam = sprintf('tp_%02i', i);
                sp = io_loadspec_sdat([mainStruct.meta.folder '\' nam '\sp\derived\' nam '_all_' sp_nam '_' condition '.SDAT'], 1);
                Cr_fwhm = op_getLW(sp, 2.8, 3.1);
                NAA_fwhm = op_getLW(sp, 1.8, 2.15);

                mainStruct.(nam).proc.(condition).(sp_nam).LWCr = Cr_fwhm;
                mainStruct.(nam).proc.(condition).(sp_nam).LWNAA = NAA_fwhm;
            end
            hp_make('save', mainStruct);
        end


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

        
        nam = sprintf('sub_%02i', id);
        % Processing steps
        % test 1: Freq and Phase aligning using ALL spectra (1) and divided
        % series (2)
        % frequency and phase aligning (NAA signal as ref)
        txt_protocol = fopen([mainStruct.meta.folder mainStruct.(nam).folder '\meta\log.txt'], 'a');
        fprintf(txt_protocol, '------ \n spectra processing steps: case - %s \n', sp_name);
        fprintf(txt_protocol, '%s \n', datetime("today"));

        pars{1, 1} = [1.9 2.1 1]; 
        sp_fpa_1 = op_freqAlignAverages_fd(sp, pars{1, 1}(1), pars{1, 1}(2), pars{1, 1}(3), 'n');
        fprintf(txt_protocol, 'op_freqAlignAverages_fd: pars - %3.1f %3.1f %3.1f \n', pars{1, 1}(:));
        sp_fpa_av = op_averaging(sp_fpa_1);  
        fprintf(txt_protocol, 'op_averaging \n');
        pars{2, 1} = {[4.4 5], 20};
        sp_fpa_av_wr = op_removeWater(sp_fpa_av, pars{2, 1}{1}, pars{2, 1}{2});
        fprintf( txt_protocol, 'op_removeWater: pars - %3.1f %3.1f %02i \n', pars{2, 1}{:});
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
        if floor(mod(mainStruct.(nam).data_check.funcTable, 10^k)/10^(k-1))>0
            ttlTime = getTTLtime([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_fmrs.xlsx']);
            regressorList = heatPain_makeRegressor([mainStruct.meta.folder mainStruct.(nam).folder '\func\' nam '_fmrs.xlsx']);
            task_starts = regressorList(:,1)-ttlTime;
            stim_dur = regressorList(:,2);
        else 
            error("There is no stimulation table data to make time_point matrix");
        end

        % using psycho py data
%         csv_file = readmatrix("G:\_other\fMRS-heatPain\_unsorted\New folder\igor_heatPainMRS_2023_Nov_10_1924.csv");
%         task_starts = csv_file(:, 19);
%         task_starts = task_starts(~isnan(task_starts));
%         task_starts = task_starts(1:37);
%         task_regr = [task_starts diff(task_starts)];
% 
%         task_starts = task_starts-task_starts(1);
        mrs_NSAmax = 315;
        mrs_timings = [1:mrs_NSAmax]*2-2;
        %need to write first active time point
        start_dynamic = 3;
        time_point_matrix = zeros(mrs_NSAmax, 1);
%         time_point_matrix(start_dynamic) = 1;
        task_starts = task_starts+4;
        k=1; series_num = 1;
        for i=start_dynamic:mrs_NSAmax
            if (task_starts(k)>mrs_timings(i)) && (task_starts(k)<mrs_timings(i)+2)
                time_point_matrix(i) = 1;
                k=k+1; series_num = 1;
                if k==length(task_starts)
                    break
                end
            else
                series_num = series_num +1;
                time_point_matrix(i) = series_num;
            end
        end
        
        mainStruct.(nam).proc.start_dynamic = start_dynamic;
        mainStruct.(nam).proc.tp_matrix = time_point_matrix;
        mainStruct.(nam).proc_check.tp_matrix = 1;
        

        %make personal HRF for heat pain stimuli during fMRS
        %simple HRF convolved with stimulus function (var task_starts)
        %get standard HRF
        MR = 16;
        dt = [0:315*MR*2]/MR;
        u = zeros(length(dt), 1);
        for i=1:length(task_starts)
            a = find(dt>task_starts(i));
            u(a(1)+ceil(dt(a(1))-task_starts(i)), 1)=1;
        end
        pars = [6, 16, 1, 1, 6, 0, 32];
        [bf, p] = spm_hrf(1/16, pars(1:6), 16);

        hrf = conv(u, bf);
        TR_samples = [0:315]*2*MR+1;
        hrf_samples=hrf(TR_samples);
        hrf_samples(1)=[];
        for i=1:max(time_point_matrix)
            hrf_mean(i) = mean(hrf_samples(time_point_matrix==i));
        end

        mainStruct.(nam).proc.hfr_mrs = hrf_mean;
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
            case 'tp_mrs'
                if mainStruct.(nam).proc_check.timepoints_spectra.sham<1
                    error('There is no time points divided data YET');
                end
                copyFiles = dir([mainStruct.meta.folder '\' nam '\sp\derived\*tp*']);
                for i = 1:length(copyFiles)
                    copyfile([copyFiles(i).folder '\' copyFiles(i).name], [out_path '\' copyFiles(i).name]);
                end

            case 'water_mrs'
                k = 2;
                if ~floor(mod(mainStruct.(nam).data_check.sp, 10^k)/10^(k-1))>0 %this takes a k-th number in data check field
                    error('there is no water reference SDAT-file or no info in mainStruct');
                end
                copyFiles = dir([mainStruct.meta.folder mainStruct.(nam).folder '\sp\' nam '_ref.*']);
                for i = 1:length(copyFiles)
                    copyfile([copyFiles(i).folder '\' copyFiles(i).name], [out_path '\' copyFiles(i).name]);
                end
                

        end

        %% Results parsing and summarising
    case 'getTableData'
        %use as hp_make('getTableData', path, id)
        mainStruct = hp_make('load');
        path = varargin{1};
        id = varargin{2};
        nam = sprintf('sub_%02i', id);
        fils = dir([path '\' nam '*\table']);
        tp_case = 0;
        for i=1:length(fils)
            temp = split(fils(i).folder, '_');
            mod_case = temp{end};
            if contains(fils(i).folder, '_tp_')
                tp_case = 1;
                tp_num = str2num(temp{end-1});
            end
            if tp_case
                out_dir = sprintf('tp_%02i_%s', tp_num, mod_case);
            else
                out_dir = sprintf('all_%s', mod_case);
            end
            mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir]);
            copyfile([fils(i).folder '\' fils(i).name], [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\' fils(i).name]);

            tp_nam = sprintf('tp_%02i', tp_num);
            mainStruct.(nam).proc.(mod_case).(tp_nam).exist = 1;
            mainStruct.(nam).proc.(mod_case).(tp_nam).path = [[mainStruct.meta.folder mainStruct.(nam).folder '\results\sp\' out_dir] '\' fils(i).name];
            mainStruct = hp_make('processLCTable',mainStruct, id, mod_case, tp_nam);
            mainStruct.(nam).proc_check.tp_spectra_res.(mod_case) = tp_num;
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

    case 'spVoxelPlacement'
        % use as hp_make('spVoxelPlacement', id)
        %makes voxel mask in MNI-space
        mainStruct = hp_make('load');
        id = varargin{1};
        nam = sprintf('sub_%02i', id);

        MRS_struct = CoRegStandAlone({[mainStruct.meta.folder mainStruct.(nam).folder '\sp\' nam '_sham.SDAT']},...
            {[mainStruct.meta.folder mainStruct.(nam).folder '\anat\' nam '_anat.nii']})
        mkdir([mainStruct.meta.folder mainStruct.(nam).folder '\anat\derived']);
        copyfile(MRS_struct.mask.vox1.outfile{1}, [mainStruct.meta.folder mainStruct.(nam).folder '\anat\derived\' nam '_sp_mask.nii']);
        delete(MRS_struct.mask.vox1.outfile{1});
        mainStruct.(nam).proc.sp_mask.path = [mainStruct.meta.folder mainStruct.(nam).folder '\anat\derived\' nam '_sp_mask.nii'];
        mainStruct.(nam).proc.sp_mask.GM = MRS_struct.out.vox1.tissue.fGM;
        mainStruct.(nam).proc.sp_mask.WM = MRS_struct.out.vox1.tissue.fWM;
        mainStruct.(nam).proc.sp_mask.CSF = MRS_struct.out.vox1.tissue.fCSF;
        mainStruct = hp_make('save', mainStruct);
        %translate to MNI
        callNormilise(mainStruct.(nam).proc.sp_mask.path,...
            [mainStruct.meta.folder mainStruct.(nam).folder '\anat\y_' nam '_anat.nii']);
        

        
    case 'getResSP'
        %use as [~, resTable] = hp_make('getResSP', condition, met)
        %gives results of the spectroscopy experiment as a matrix
        mainStruct = hp_make('load');
        
        condition = varargin{1};
        met = varargin{2};
        
        if contains(met, 'LW')
            resTable.LW.Cr = zeros([mainStruct.meta.subNumbers, 6]);
            resTable.LW.NAA = zeros([mainStruct.meta.subNumbers, 6]);
        else
            resTable.(met).Conc = zeros([mainStruct.meta.subNumbers, 6]);
            resTable.(met).ConcCr = zeros([mainStruct.meta.subNumbers, 6]);
        end

        for id=1:mainStruct.meta.subNumbers
            nam = sprintf('sub_%02i', id);
            if mainStruct.(nam).proc_check.tp_spectra_res.(condition)>0
                colNum = mainStruct.(nam).proc_check.tp_spectra_res.(condition);
                for ii=1:colNum
                    tp_nam = sprintf('tp_%02i', ii);
                    if contains(met, 'LW')
                        resTable.LW.NAA(id, ii) = mainStruct.(nam).proc.(condition).(tp_nam).LWNAA;
                        resTable.LW.Cr(id, ii) = mainStruct.(nam).proc.(condition).(tp_nam).LWCr;
                    else
                        resTable.(met).Conc(id, ii) = mainStruct.(nam).proc.(condition).(tp_nam).Glx;
                    end

                end
            end
        end
        varargout{1} = resTable;
                

%% Q6_time_series_task

    case 'Q6'

        % to activate the case -    hp_make('Q6',8)
        
        id = varargin{1};
        mainStruct = hp_make('load');
        nam = sprintf('sub_%02i', id);

        fmri_file = ([mainStruct.meta.folder '\' nam '\func\' nam '_func.nii']);
        rp_file = ([mainStruct.meta.folder '\' nam '\derived\rp_' nam '_func.txt']);
        SPM_file = ([mainStruct.meta.folder '\' nam '\derived\res\SPM.mat']);
        
        %fmri_file = ('C:\Helen\Docs\DHiT\Alex\sub_08\func\sub_08_func.nii');
        %rp_file = ('C:\Helen\Docs\DHiT\Alex\sub_08\derived\rp_sub_08_func.txt');
        %SPM_file = ('C:\Helen\Docs\DHiT\Alex\sub_08\derived\res\SPM.mat');

        Q6_time_series_task(fmri_file, rp_file, SPM_file);
        saveas(gcf,[mainStruct.meta.folder '\' nam '\meta\Q6.jpg'],'jpg');
         



        %% Some not very important features or minor experiments (consider to remove it)
        
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

end
end

function startTime = getTTLtime(xlsfile)
% stimulus times ONLY for MRS modality
    eventsTable = readtable(xlsfile);
    TTL = strfind(eventsTable.Events, 'TTL');
    for i=1:length(TTL)
        if ~isempty(TTL{i})
            startTime = eventsTable.Timestamp_msec_(i)/1000;
            break
        end
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