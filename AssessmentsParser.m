function [rating, reaction_time] = AssessmentsParser(filename)
% read log-file
log_table = fopen(filename, 'r');
i=1;
s1{i} = fgetl(log_table);
while s1{i}>0
    i=i+1;
    s1{i} = fgetl(log_table);
end

k_exp = 1;
k_ans = 1;
for i=1:length(s1)
    if isstr(s1{i})
        if contains(s1{i}, 'EXP 	New trial')
            temp = split(s1{i}, ' ');
            exp_time{k_exp, 1} = str2num(temp{1});
            exp_time{k_exp, 2} = sscanf(temp{5}, '(rep=%i,');
            k_exp = k_exp + 1;
        end
        if contains(s1{i}, 'DATA 	Keypress: b')
            temp = split(s1{i}, ' ');
            ans_time{k_ans, 1} = str2num(temp{1});
            k_ans = k_ans+1;
        end
        if contains(s1{i}, 'DATA 	Keypress: a')
            temp = split(s1{i}, ' ');
            ans_time{k_ans, 2} = str2num(temp{1});
            k_ans = k_ans+1;
        end
    end
end


for i=1:length(ans_time)
    if ~isempty(ans_time{i, 1})
        key_resp(i, 1) = ans_time{i, 1};
    end
    if ~isempty(ans_time{i, 2})
        key_resp(i, 1) = ans_time{i, 2};
    end
end

k=1;
reaction_time = zeros(50, 1);
for i=1:length(exp_time)
    temp_resp = key_resp - exp_time{i, 1};
    reaction_time(i) = min(temp_resp(temp_resp>0));
end

reaction_time(reaction_time>10) = [];

b_count = zeros(50, 1)-1;
a_count = zeros(50, 1);
rating = zeros(50,1)+4;
for i=1:length(exp_time)-1
    for j=1:length(ans_time)
        if (key_resp(j)>exp_time{i, 1}) && (key_resp(j)<exp_time{i+1, 1})
            if ~isempty(ans_time{j, 1})
                b_count(i)=b_count(i)+1;
            else
                a_count(i)=a_count(i)+1;
            end
        end
    end
    rating(i) = rating(i)+b_count(i)-a_count(i);
end
end