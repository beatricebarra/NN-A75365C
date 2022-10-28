%%
% @Author: barrab
% @Date:   2019-06-13 12:14:23

%% Loading dataset 


% Load data
datapath = 'datapath'; 
Elbow=  readtable(fullfile(datapath, 'Figure3B_ElbowExtension.csv')); 
Hand=  readtable(fullfile(datapath, 'Figure3B_HandFlexion.csv')); 


% Load colors
OneTrueCM; 


%% Plot Figure 3B

figure
frequencies = Elbow.Properties.VariableNames'; 

nFreq = length(frequencies); 
subplot(1, 2, 1)
hold on
for iF = 2 : nFreq
    datapoints = Elbow.(frequencies{iF}); 
    
    plot(iF*ones(size(datapoints))+ 0.2*(rand(length(datapoints),1)-0.5), datapoints, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, nanmedian(datapoints), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

ylim([0 60])
ylabel('degrees')
set(gca, 'xtick', [1 : length(frequencies)], 'xticklabel', frequencies)
title('Elbow excursion')

subplot(1, 2, 2)
frequencies = Hand.Properties.VariableNames'; 
hold on
for iF = 2 : nFreq
    datapoints = Hand.(frequencies{iF}); 
    
    plot(iF*ones(size(datapoints))+ 0.2*(rand(length(datapoints),1)-0.5), datapoints, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
     plot(iF, nanmedian(datapoints), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

ylabel('degrees')
set(gca, 'xtick', [1 : length(frequencies)], 'xticklabel', frequencies)
title('Hand Flexion')