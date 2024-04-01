function [rating, reaction_time] = AssessmentParser2(filename)
% reads csv. file from psychopy and gives pain response timing and pain
% assessments
 ppy_table = readtable(filename);

 %get reaction time
 reaction_time = ppy_table.key_resp_2_rt';
 reaction_time(isnan(reaction_time)) = [];

%get pain assessments
rating = []; k =1;
rat_hist = ppy_table.rating_history;
for i=1:length(rat_hist)
    if isempty(rat_hist{i}) || contains(rat_hist{i}, 'None')
        continue
    end
    temp = split(rat_hist{i}, ',');
    temp{end-1} = regexp(temp{end-1}, '\d', 'match');
    rating(k, 1) = str2num(temp{end-1}{1});
    k = k+1;
end
rating = rating(5:5:end);