warning('off','all')
clear all
clc
close all
Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES/';
load('/Users/barrab01/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')


Animal = 'Brienne';
expDates = { '20190529', '20190611', '20190613', '20190618', '20190624'};

selectedFiles = {
    [5:8 14:17]; ...
    [1,  8, 12];...
    [2 11];... 8 %[1 3 4 5 13];... 8 
    [2, 3,];
    [1 2  5  8  13  18];...
}
stimMask = {
    [2 2 2 2 2 2 2 2 ]; ...
    [0, 0,0]; ... 
    [0 0]; ...
    [0 0 ];
    [0 0 0 0 0 0]; ... 
}% each cell is a date, 1 means stim file, 0 is basel

scspapercolors = [230 230 230; 125 125 125; 255 183 27];
recsystems = {'TDT', 'VICON', 'BLACKROCK'};
additionalpath = 'TDT';

set(0, 'DefaultFigureRenderer', 'painters');
iif = 0;

% Retrieve intact baselines count from files

iDate = 1;

intact_data = [];
for iF = 1 : length(selectedFiles{iDate})
   
    expDate = expDates{iDate};
    filenameVICON = [expDate,'_',Animal, '_'];
    dataPath = fullfile(Root, Animal,expDate, additionalpath);

    iFile = selectedFiles{iDate}(iF);
    disp(['Loading File n ', num2str(iFile) ])

    if any(strcmp(recsystems,'VICON'))
        VICONfile = [filenameVICON, num2str(iFile, '%02.f')];
        VICONdat{iDate}{iF} = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '.mat' ]));
    end
    
    Good_trials = find(cellfun(@(x) strcmp(x, 'EndTrial'), VICONdat{iDate}{iF} .event.labels));
    nGoodtrials{iDate}(iF) = length(Good_trials);
    time_trial{iDate}(iF) = size( VICONdat{iDate}{iF}.analog.data, 1)*(1/1000);
    
    goodtrials_vect = zeros(1, ceil(time_trial{iDate}(iF))); 
    goodtrials_vect(1:nGoodtrials{iDate}(iF)) = 1; 
    intact_data = [intact_data, goodtrials_vect];
    
end

nSamples = 1000;
   
[mean_event_I, std_event_I] =  compute_stats_bootstrap (intact_data, nSamples); 
   
    
%% Retrieve lesion data

Animal = 'Brienne';
expDates = {  '20190611', '20190613', '20190618', '20190624'};

selectedFiles = {
    [1,  8, 12];...
    [2 11];... 8 %[1 3 4 5 13];... 8 
    [2, 3,];
    [1 2  5  8  13  18];...
}
stimMask = {
   
    [0, 0,0]; ... 
    [0 0]; ...
    [0 0 ];
    [0 0 0 0 0 0]; ... 
}% each cell is a date, 1 means stim file, 0 is basel


timefile = xlsread('/Users/barrab01/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/REVIEW-ROUND2/Events_data/ActiveTimesMkBr.xlsx'); 
dates_file = cellstr(num2str(timefile(:,1))); 
trials_file = timefile(:,2); 

active_times = cell(length(expDates ), 1);
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

for iD = 1 : length(events)
    baseline_data.reach = [];
    baseline_data.grasp = [];
    baseline_data.pull = [];

    for iT = 1 : length(events{iD})
        times = active_times{iD}{iT}; 
        active_bins = times.stops-times.starts; 
        trial_time{iD}(iT)= sum(active_bins(not(isnan(active_bins))));
        
        Reach_num{iD}(iT) = events{iD}{iT}.reach; 
        Grasp_num{iD}(iT) = events{iD}{iT}.grasp; 
        Pull_num{iD}(iT) = events{iD}{iT}.pull;
        
        reach_vect = zeros(1, ceil(trial_time{iD}(iT))); 
        try reach_vect(1:Reach_num{iD}(iT)) = 1; catch end
        grasp_vect = zeros(1, ceil(trial_time{iD}(iT))); 
        try grasp_vect(1:Grasp_num{iD}(iT)) = 1;catch end
        pull_vect = zeros(1, ceil(trial_time{iD}(iT))); 
        pull_vect(1:Pull_num{iD}(iT)) = 1; 
        
        
        baseline_data.reach = [baseline_data.reach, reach_vect];
        baseline_data.grasp = [baseline_data.grasp,grasp_vect ];
        baseline_data.pull = [baseline_data.pull,pull_vect ];
    end
    
    
    % Analysis of rate
    nSamples = 1000;
   
    [mean_reach_B(iD), std_reach_B(iD)] =  compute_stats_bootstrap (baseline_data.reach, nSamples); 
    
    [mean_grasp_B(iD), std_grasp_B(iD)] =  compute_stats_bootstrap (baseline_data.grasp, nSamples); 
    
    [mean_pull_B(iD), std_pull_B(iD)] =  compute_stats_bootstrap (baseline_data.pull, nSamples); 

    disp ('    ') 
    
end


%%
date_idxs = [1:length(expDates )];
figure
subplot(1, 3, 1) % Reach ration
hold on
plot(0, mean_event_I, '-ob')
plot([0 0], [mean_event_I-std_event_I, mean_event_I+std_event_I], 'b')

plot_meanswithstds(date_idxs, mean_reach_B, std_reach_B, 'k')

%set(gca, 'xtick', date_idxs, 'xticklabels', expDates)
%xlim([date_idxs(1) - 1, date_idxs(end)+1])
title('Reach rate')

subplot(1, 3, 2) % Reach ration
hold on
plot(0, mean_event_I, '-ob')
plot([0 0], [mean_event_I-std_event_I, mean_event_I+std_event_I], 'b')

plot_meanswithstds(date_idxs, mean_grasp_B, std_grasp_B, 'k')
title('Grasp rate')


subplot(1, 3, 3) % Reach ration
hold on
plot(0, mean_event_I, '-ob')
plot([0 0], [mean_event_I-std_event_I, mean_event_I+std_event_I], 'b')

plot_meanswithstds(date_idxs, mean_pull_B, std_pull_B, 'k')

