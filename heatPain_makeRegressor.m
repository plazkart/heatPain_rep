%
function  regressorList = heatPain_makeRegressor(taskFile)
 xlFile = xlsread(taskFile);

time = xlFile(:, 1);
temp = xlFile(:, 3);

tempStart = min(temp(1:1000))+ (max(temp) - min(temp(1:1000)))/2;

k=1;
temp = abs(temp - tempStart);

temp_minimals = temp<0.8;
idxs = 1:length(temp);
idxs = idxs(temp_minimals);

startTimePoints = [];


while length(idxs)>1
    idxs_local = idxs((time(idxs)-time(idxs(1)))<2000);
    [~, idxMin] = min(temp(idxs_local));
    startTimePoints(k) = time(idxs(1)+idxMin-1);
    k=k+1;
    idxs = setdiff(idxs, idxs_local);
end

startTimePoints = startTimePoints/1000;
if mod(length(startTimePoints), 2)>0
    startTimePoints(end) = [];
end
regressorList = [startTimePoints(1:2:end)', startTimePoints(2:2:end)'-startTimePoints(1:2:end)' ];

%% baseline regressor list
if false
    % rgList2 = [0 startTimePoints(1); reshape([startTimePoints(2:2:52) startTimePoints(3:2:end)-startTimePoints(2:2:52)], 26, 2); startTimePoints(54) 12];
    cellOut = {'onset', 'duration'};
    for i=1:size(regressorList, 1)
        cellOut{i+1, 1} = regressorList(i,1);cellOut{i+1, 2} = regressorList(i, 2);
    end
    % writematrix(regressorList,'G:\_other\fMRI-thermal\_bids\sub_010\ses_02\func\sub-010_ses_02_task-hp_events.tsv', 'filetype','text', 'Delimiter','tab');
    writecell(cellOut,'G:\_other\fMRI-thermal\_bids\sub-010\ses-02\func\sub-010_ses-02_task-hp_events.tsv', 'filetype','text', 'Delimiter','tab');
end