warning('off','all')
clear all
clc
close all
% Root = 'P:\data_raw\primate\SCS_Monkey\';
% currentpath = 'C:\GIT\unifr_cervicalEES\';
% load('C:\GIT\unifr_cervicalEES\Structs\onetrueCM.mat')


currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES/';
load('/Users/barrab01/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')
% 10 06
% [2 3 10 11 12  15 16 17 20]; ...
% [0 0 1 1 1 1 0 1 1]; ...
Animal = 'Brienne';
expDates = {'20190611', '20190613', '20190618', '20190624'};
selectedFiles = {
    [1, 5, 7, 8, 10, 12];...
    [2 6 7 11];... 8 %[1 3 4 5 13];... 8 
    [2, 3, 5, 6, 7];
    [1 2 3 4 5 6 7 8 10 11 12 13 14 15 18];...
}
stimMask = {
    [0, 1, 1, 0, 1 ,0]; ... 
    [0 1 1 0]; ...
    [0 0 1 1 1];
    [0 0 1 1 0 1 1 0 1 1 1 0 1 1 0]; ... 
}% each cell is a date, 1 means stim file, 0 is basel



timefile = xlsread('/Users/barrab01/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/REVIEW-ROUND2/Events_data/ActiveTimesMkBr.xlsx'); 
dates_file = cellstr(num2str(timefile(:,1))); 
trials_file = timefile(:,2); 

active_times = cell(length(expDates), 1);
for iD = 1 : length(expDates)
    active_times{iD} = cell(length(selectedFiles{iD}), 1);
    for iT = 1 : length(active_times{iD})
        active_times{iD}{iT}.starts = [];
        active_times{iD}{iT}.stops = [];
    end
end

for i = 1 : size(timefile, 1)
    try
        idx_date = find(cellfun(@(x) strcmp(x, dates_file{i}), expDates));
    
        idx_trial = find(selectedFiles{idx_date} == trials_file(i)); 
    
        last_col = size(timefile, 2); 
        total_times{idx_date}(idx_trial)= timefile(i,last_col);
        for iC = 3 : 2 : size(timefile, 2)-1 
            if not(isnan(timefile(iC)))
                active_times{idx_date}{idx_trial}.starts = ...
                    [active_times{idx_date}{idx_trial}.starts, timefile(i, iC)];
                active_times{idx_date}{idx_trial}.stops = ...
                    [active_times{idx_date}{idx_trial}.stops, timefile(i, iC + 1)];
            end
        end
    catch
    end
            
end



events_file = xlsread('/Users/barrab01/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/REVIEW-ROUND2/Events_data/MkBr_events.xlsx'); 
dates_file_ev = cellstr(num2str(events_file(:,1))); 
trials_file_ev = events_file(:,2); 
events = cell(length(expDates), 1);
for iD = 1 : length(expDates)
    events{iD} = cell(length(selectedFiles{iD}), 1);
    for iT = 1 : length(active_times{iD})
        events{iD}{iT}.reach = [];
        events{iD}{iT}.grasp = [];
        events{iD}{iT}.pull = [];
    end
end

for i = 1 : size(events_file, 1)
    try
        idx_date = find(cellfun(@(x) strcmp(x, dates_file_ev{i}), expDates));
        idx_trial = find(selectedFiles{idx_date} == trials_file_ev(i)); 
    
        events{idx_date}{idx_trial}.reach = ...
            [events{idx_date}{idx_trial}.reach, events_file(i, 3)];
        events{idx_date}{idx_trial}.grasp = ...
            [events{idx_date}{idx_trial}.grasp, events_file(i, 4)];
        
        events{idx_date}{idx_trial}.pull = ...
            [events{idx_date}{idx_trial}.pull, events_file(i, 5)];    
    
    catch
    end
            
end


baseline_data.reach = [];
baseline_data.grasp = [];
baseline_data.pull = [];

stim_data.reach = [];
stim_data.grasp = [];
stim_data.pull = [];


for iD = 1 : length(events)
    baseline_data.reach = [];
    baseline_data.grasp = [];
    baseline_data.pull = [];

    stim_data.reach = [];
    stim_data.grasp = [];
    stim_data.pull = [];
    for iT = 1 : length(events{iD})
        times = active_times{iD}{iT}; 
        active_bins = times.stops-times.starts; 
        trial_time{iD}(iT)= sum(active_bins(not(isnan(active_bins))));
        
        Reach_num{iD}(iT) = events{iD}{iT}.reach; 
        Grasp_num{iD}(iT) = events{iD}{iT}.grasp; 
        Pull_num{iD}(iT) = events{iD}{iT}.pull;
        
       
        
        
        Reach_rate{iD}(iT) = events{iD}{iT}.reach/trial_time{iD}(iT); 
        Grasp_rate{iD}(iT) = events{iD}{iT}.grasp/trial_time{iD}(iT); 
        Pull_rate{iD}(iT) = events{iD}{iT}.pull/trial_time{iD}(iT); 
        
        reach_vect = zeros(1, ceil(trial_time{iD}(iT))); 
        try reach_vect(1:Reach_num{iD}(iT)) = 1; catch end
        grasp_vect = zeros(1, ceil(trial_time{iD}(iT))); 
        try grasp_vect(1:Grasp_num{iD}(iT)) = 1;catch end
        pull_vect = zeros(1, ceil(trial_time{iD}(iT))); 
        pull_vect(1:Pull_num{iD}(iT)) = 1; 
        iD, iT
        find(pull_vect)
        
        if stimMask{iD}(iT)==0
            baseline_data.reach = [baseline_data.reach, reach_vect];
            baseline_data.grasp = [baseline_data.grasp,grasp_vect ];
            baseline_data.pull = [baseline_data.pull,pull_vect ];
        elseif stimMask{iD}(iT)==1
            stim_data.reach = [stim_data.reach, reach_vect];
            stim_data.grasp = [stim_data.grasp,grasp_vect ];
            stim_data.pull = [stim_data.pull,pull_vect ];
        end
    end
    
    
    % Analysis of rate
    nSamples = 1000;
   
    [mean_reach_B(iD), std_reach_B(iD)] =  compute_stats_bootstrap (baseline_data.reach, nSamples); 
    [mean_reach_S(iD), std_reach_S(iD)] =  compute_stats_bootstrap (stim_data.reach, nSamples); 
    
    [mean_grasp_B(iD), std_grasp_B(iD)] =  compute_stats_bootstrap (baseline_data.grasp, nSamples); 
    [mean_grasp_S(iD), std_grasp_S(iD)] =  compute_stats_bootstrap (stim_data.grasp, nSamples); 
    
    [mean_pull_B(iD), std_pull_B(iD)] =  compute_stats_bootstrap (baseline_data.pull, nSamples); 
    [mean_pull_S(iD), std_pull_S(iD)] =  compute_stats_bootstrap (stim_data.pull, nSamples); 
    
    
    
    
    disp ('    ') 
    
end


%
date_idxs = [1 2 3 4];
figure
subplot(2, 3, 1) % Reach ration
hold on
plot_meanswithstds(date_idxs, mean_reach_B, std_reach_B, 'k')
plot_meanswithstds(date_idxs, mean_reach_S, std_reach_S, 'r')

for i = 1 : length(date_idxs)
    [p] = manual_pvalue(mean_reach_S(i), mean_reach_B(i), std_reach_S(i) , std_reach_B(i))
    if p <= 0.05
        text(i, (mean_reach_S(i) + mean_reach_B(i))/2, '*')
        %set(gca,'fontsize', 18)
    end
end

set(gca, 'xtick', date_idxs, 'xticklabels', expDates)
xlim([date_idxs(1) - 1, date_idxs(end)+1])
title('Reach /attempts in different dates')

subplot(2, 3, 2) % Grasp
hold on
plot_meanswithstds(date_idxs, mean_grasp_B, std_grasp_B, 'k')
plot_meanswithstds(date_idxs, mean_grasp_S, std_grasp_S, 'r')

for i = 1 : length(date_idxs)
    [p] = manual_pvalue(mean_grasp_S(i), mean_grasp_B(i), std_grasp_S(i) , std_grasp_B(i))
    if p <= 0.05
        text(i, (mean_grasp_S(i) + mean_grasp_B(i))/2, '*')
        %set(gca,'fontsize', 18)
    end
end
title('Grasp /attempts in different dates')
set(gca, 'xtick', date_idxs, 'xticklabels', expDates)
xlim([date_idxs(1) - 1, date_idxs(end)+1])

subplot(2, 3, 3) % Pull
hold on
plot_meanswithstds(date_idxs, mean_pull_B, std_pull_B, 'k')
plot_meanswithstds(date_idxs, mean_pull_S, std_pull_S, 'r')

for i = 1 : length(date_idxs)
    
    [p] = manual_pvalue(mean_pull_S(i), mean_pull_B(i), std_pull_S(i) , std_pull_B(i))
    if p <= 0.05
        text(i, (mean_pull_S(i) + mean_pull_B(i))/2, '*')
        %set(gca,'fontsize', 18)
    end
end
title('Pull / attemptsnumber in different dates')
set(gca, 'xtick', date_idxs, 'xticklabels', expDates)
xlim([date_idxs(1) - 1, date_idxs(end)+1])


%%
%%
tot_b_reaches = 0 ; 
tot_b_grasps= 0 ; 
tot_b_pulls= 0 ; 
tot_s_reaches = 0 ; 
tot_s_grasps= 0 ; 
tot_s_pulls= 0 ; 
for iD = 1 : length(events)
    for iT = 1 : length(events{iD})
       
        if stimMask{iD}(iT)==0
            tot_b_reaches = tot_b_reaches + Reach_num{iD}(iT);
            tot_b_grasps = tot_b_grasps + Grasp_num{iD}(iT);
            tot_b_pulls = tot_b_pulls + Pull_num{iD}(iT);
        else
            tot_s_reaches = tot_s_reaches + Reach_num{iD}(iT);
            tot_s_grasps = tot_s_grasps+ Grasp_num{iD}(iT);
            tot_s_pulls = tot_s_pulls + Pull_num{iD}(iT); 
        end
        
    end
end


figure
subplot(1, 3, 1)
pie([tot_b_reaches, tot_s_reaches])
subplot(1, 3, 2)
pie([tot_b_grasps, tot_s_grasps])
subplot(1, 3, 3)
pie([tot_b_pulls, tot_s_pulls])