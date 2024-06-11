function Q6_time_series_task(fmri_file, rp_file, SPM_file, ses_num)
% clear
% 
% subjfolder = dir('F:\fMRI_QC\fmri-open-qc-task_bias_corrected\sub*');

% figure('WindowState','maximized');
figure('position',[10 50 1500 500]);
% for subji = 1:length(subjfolder)
    subji = 1;
    fprintf('Processing subject %d ...',subji)
    
    % load the raw functional images and calculate the global mean and pairwise variance
    clear a dt
    img_number = spm_vol(fmri_file);
    img_number = length(img_number);
    
    for imagei = 1:img_number
%         v = spm_vol(fullfile('F:\fMRI_QC\fmri-open-qc-task_bias_corrected\', subjfolder(subji).name, 'func', [subjfolder(subji).name '_task-pamenc_bold.nii,' num2str(imagei)]));
v = spm_vol([fmri_file ',' num2str(imagei)]);
        y = spm_read_vols(v);
        
        a(:,imagei) = y(:);
    end
    
    gm = mean(mean(a));
    
    % calculate pairwise variance
    for imagei = 1:img_number-1
        dt(imagei) = (mean((a(:,imagei) - a(:,imagei+1)).^2))/gm;
    end
    
    meany = mean(a)./gm;
    
    
    % load rigid body motion parameters and calculate framewise displacement 
    rp = load(rp_file);

    fd_trans = fd_calc(rp(:,1:3));
    fd_rotat = fd_calc(rp(:,4:6)*180/pi);
    
    fd_max_trans(subji,1) = max(fd_trans);
    fd_max_rotat(subji,1) = max(fd_rotat);
    
    fd_mean_trans(subji,1) = mean(fd_trans);
    fd_mean_rotat(subji,1) = mean(fd_rotat);
    
    
    % load task design
    load(SPM_file);
  
    
    % plots
    subplot(2,3,1)
    plot(meany)
    title('Global mean (raw)');xlabel('Image number')
    box off
    
    subplot(2,3,4)
    plot(dt)
    yline(mean(dt)+3*std(dt),'-','3 SD','color',[0 0.4470 0.7410]);
    title('Pairwise variance (raw)');xlabel('Image pair')
    box off
    
    subplot(2,3,2)
    plot([rp(:,1:3) rp(:,4:6)*180/pi])
    title('Rigid body motion');xlabel('Image number')
    box off
    
    subplot(2,3,5)
    plot([fd_trans fd_rotat])
    title('Framewise displacement');xlabel('Image pair')
    legend('Translation','Rotation','location','best','box','off')
    if max(max([fd_trans fd_rotat])) > 1.5
        yline(1.5,'r','FD = 1.5');
    end
    box off
    
%     subplot(2,4,3)
%     plot(SPM.xX.X(:,[1, end]))
%     title('Task design');xlabel('Image number')
%     legend('Task','Control','location','best','box','off')
%     box off
%     
%     subplot(2,4,7)
%     plot(diff(SPM.xX.X(:,[1, end])))
%     title('d(Task)');xlabel('Image pair')
%     legend('Task','Control','location','best','box','off')
%     box off
    
    sess_number = size(SPM.Sess, 2);
    dm_size = size(SPM.xX.X);
    if nargin<4
        ses_num = 1;
    end
    idxs = [(ses_num-1)*dm_size(1)/sess_number+1 : (ses_num)*dm_size(1)/sess_number];
    idxs_col = [(ses_num-1)*dm_size(2)/sess_number+1 : (ses_num)*dm_size(2)/sess_number];
    
    design_matrix = SPM.xX.X(idxs, idxs_col);
    corr_raw = corr([meany' rp design_matrix(:,1)]);
    corr_dt = corr([dt' fd_trans fd_rotat diff(design_matrix(:,1))]);
    
    subplot(2,3,3)
    imagesc(corr_raw)
    caxis([-1 1]);
    yticks(1:9)
    yticklabels({'G raw','HM 1','HM 2','HM 3','HM 4','HM 5','HM 6','Task','Control'})
    colorbar
    title('Correlation');
    %mat_mask1 = [zeros(7,7) ones(7,1);ones(1,7) zeros(1,1)];
    mat_mask1 = [zeros(1,1) ones(1,7);zeros(7,7) zeros(7,1)];
    disp(mat_mask1);
    [x,y] = find((abs(corr_raw)>0.7).*mat_mask1);
    hold on; scatter(x,y,[],'r','filled')
    
    subplot(2,3,6)
    imagesc(corr_dt)
    caxis([-1 1]);
    yticks(1:5)
    yticklabels({'V raw','FD T','FD R','d(Task)','d(Control)'})
    colorbar
    title('Correlation');
    mat_mask2 = [0 1 1 1 ;0 0 1 1 ;0 0 0 1 ;0 0 0 0 ];
    disp(mat_mask2);
    [x,y] = find((abs(corr_dt)>0.7).*mat_mask2);
    hold on; scatter(x,y,[],'r','filled')

%     exportgraphics(gcf,['C:\Users\synge\Documents\GitHub\Preprocessing_and_QC\QC_images\Q6_time_series_task\' subjfolder(subji).name '.jpg']);
    
    fprintf('done!\n')
end


function fd = fd_calc(Y)
% calculate framewise displacement
    Y_diff = diff(Y);
    multp = Y_diff*Y_diff';
    fd = sqrt(diag(multp));
end