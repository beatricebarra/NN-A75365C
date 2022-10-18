%%
% @Author: barrab
% @Date:   2019-06-13 12:14:23

%% Loading dataset 


% Load data
datapath = '/Users/barrab01/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/EDITOR ROUND/SourceData/source_data_final_ODC/'; 
E1=  readtable(fullfile(datapath, 'Figure2A_E1_Mk-Yg.csv')); 
E4=  readtable(fullfile(datapath, 'Figure2A_E4_Mk-Yg.csv')); 
E6=  readtable(fullfile(datapath, 'Figure2A_E6_Mk-Yg.csv')); 

% Load colors
OneTrueCM; 

%% Plotting Recruitment curves

h =figure('units','normalized','outerposition',[0 0 0.5 1]);
% Plot Rostral site
nChan = size(E1, 2) -1;
labels = E1.Properties.VariableNames'; 
muscles = labels(2:end, 1); 
subplot(3, 1, 1)
hold on
for ch = 1 : nChan
    x = E1.Current; 
    y = E1.(labels{ch+1}); 
    plot(x, y , 'Linewidth', 2, 'Color', muscleColors.(labels{ch+1}))
end
ylim([0 1])
title('E1')
legend(muscles ,'Location','NorthWest')

% Plot Medial site
nChan = size(E4, 2) -1;
labels = E4.Properties.VariableNames'; 
muscles = labels(2:end, 1); 
subplot(3, 1, 2)
hold on
for ch = 1 : nChan
    x = E4.Current; 
    y = E4.(labels{ch+1}); 
    plot(x, y , 'Linewidth', 2, 'Color', muscleColors.(labels{ch+1}))
end
ylim([0 1])
title('E1')
legend(muscles ,'Location','NorthWest')

% Plot Caudal site
nChan = size(E6, 2) -1;
labels = E6.Properties.VariableNames'; 
muscles = labels(2:end, 1); 
subplot(3, 1, 3)
hold on
for ch = 1 : nChan
    x = E6.Current; 
    y = E6.(labels{ch+1}); 
    plot(x, y , 'Linewidth', 2, 'Color', muscleColors.(labels{ch+1}))
end
ylim([0 1])
title('E1')
legend(muscles ,'Location','NorthWest')
