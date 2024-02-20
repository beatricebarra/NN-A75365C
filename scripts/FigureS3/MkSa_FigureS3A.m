%%
% @Author: barrab
% @Date:   2019-06-13 12:14:23
% @Last Modified by:   barrab
% @Last Modified time: 2019-06-13 12:14:23

% options.format = 'pdf'; options.showCode = false;publish('Ygritte',options);
%% Recruitment curves have been performed on : 
% 
% * 06.08.2019
% * 14.08.2019
% * 19.08.2019
% * 26.08.2019

clear all
close all


pittPC =1;
Animal = 'Sansa';

params_18092018;


nResp = 4;
nChan = 3;
fS = 12207.03;
lengthResp = floor(0.05*fS);
singleLead = 0;
OneTrueCM;
emnmappingsansa;
%load('/Users/Bea/Documents/PhD/Repos/unifr_cervicalEES/Structs/papercolors.mat')
load('C:\GIT\unifr_cervicalEES\Structs\papercolors.mat')
%% Loading files ...

for i = 1 : length(fileNames)
    
    % Loading
    TDTdat{i} = load(fullfile(Root, [fileNames{i}, '.mat']));
    % Saving additional infos
    TDTdat{i}.currents = currents{i};
    TDTdat{i}.muscles = muscles;
end  
%% .. and preprocessing
for i = 1 : length(fileNames)
    % Deriving bipolar EMGs from monopolar acquisitions
    TDTdat{i}.EMGtb = zeros(size(TDTdat{i}.EMGt,1), nChan);
    if singleLead
       for ch = 1:nChan
            TDTdat{i}.EMGtb(:,ch) = TDTdat{i}.EMGt(:,emgmapping(ch, 2));
        end 
    else
        for ch = 1:1:(nChan)
            TDTdat{i}.EMGtb(:,ch) = TDTdat{i}.EMGt(:,emgmapping(ch, 1))-TDTdat{i}.EMGt(:,emgmapping(ch, 2));
        end
    end
    if i == 6 
        TDTdat{i}.EMGtb = [TDTdat{i}.EMGtb(1:610, :); TDTdat{i}.EMGtb];
    end
        
    % Number of explored currents during recruitment curves
    %nCurrents = size(TDTdat{i}.EMGtb, 1)/(lengthResp * nResp);
    nCurrents = length(TDTdat{i}.currents);
    % Computing recruitment curves based on peak to peak measure
    for ch = 1 : nChan
        for iCurr = 1 : nCurrents
            allResps = reshape(TDTdat{i}.EMGtb((iCurr-1)*lengthResp*nResp+1: (iCurr)*lengthResp*nResp, ch), [lengthResp, nResp]);
            allResps(1:50,:) = 0;
            subplot(6, 3, 3*(i-1) + ch)
            plot(allResps)
            
            avgResp = median(allResps, 2);
            TDTdat{i}.trigResp{ch}(:,:,iCurr) = allResps;
            TDTdat{i}.p2p(ch, iCurr) = max(avgResp) - min(avgResp);
        end
    end
    pause
    
    
    
end

% Normalizing recruitment curves on all the trials    
maxResp = zeros(nChan,1);
for ch = 1 : nChan
    for iF = 1 : length(fileNames)
        maxResp(ch, 1) = max([maxResp(ch, 1), max(TDTdat{iF}.p2p(ch,:))]);
    end
end

for iF = 1 : length(fileNames)  
    for ch = 1 : nChan
        TDTdat{iF}.RC(ch,:) = TDTdat{iF}.p2p(ch,:)./maxResp(ch, 1) ;
        %TDTdat{iF}.stdRC(ch,:) = std(TDTdat{i}.trigResp')';
    end
end
maxactivation = 0;
minactivation  = 1000;
for iM = 1 : length(TDTdat{iF}.muscles)
    for iF = 1 : length(fileNames)       
        maxactivation(iM) = max([max(TDTdat{iF}.trigResp{iM}(:)), maxactivation]);
        minactivation(iM) = min([min(TDTdat{iF}.trigResp{iM}(:)), minactivation]);
    end
end

MAXACT = max(maxactivation);
MINACT = min(minactivation);
snipLength_time = (1/fS)*lengthResp;
dt = (1/fS);
dims = numSubplots(length(pinNames));
maxactivation = 0 ; 
minactivation = 1000;
for iM = 1 : length(TDTdat{iF}.muscles)  
    time(iM, :) = linspace((iM-1)*snipLength_time, iM*(snipLength_time)-dt, lengthResp);
end

for iF = 1 : length(fileNames)
    RCdat{iF}.currents = currents{iF};
    RCdat{iF}.electrode = pinNames{iF};
    RCdat{iF}.muscles = TDTdat{iF}.muscles;
    RCdat{iF}.data = TDTdat{iF}.RC;

end
%save('C:\GIT\unifr_cervicalEES\Structs\RCSansa', 'RCdat')





%% Plotting Recruitment curves

h =figure('units','normalized','outerposition',[0 0 0.5 1]);
dims = numSubplots(length(fileNames) );
for iF = 1 : length(fileNames)  
    subplot(dims(1), dims(2), iF)
    hold on
    for ch = 1 : nChan
        plot(currents{iF}, TDTdat{iF}.RC(ch,:), 'Linewidth', 2, 'Color', muscleColors.(muscles{ch}))
    end
    ylim([0 1])
    title(pinNames{iF})
    legend(muscles ,'Location','NorthWest')
end


%% Plotting activation index as the motor pools


musclexposition = [ 2 4 7 ]/4;
    %]/2;
electrodeyposition = [6 5 4 3 2 1 ];
h =figure('units','normalized','outerposition',[0 0 0.5 1]);
hold on
clear activationindex

for ch = 1 : 3
    for iF = 1 : length(fileNames) 
        activationindex(iF) = mean(TDTdat{iF}.RC(ch,:));
        plot(activationindex(iF) + musclexposition(ch), electrodeyposition(iF),  'o', 'Linewidth', 2,...
            'MarkerFaceColor', muscleColors.(muscles{ch}), ...
            'MarkerEdgeColor', muscleColors.(muscles{ch}), ...
            'Markersize', 10)
    end
    
    plot(activationindex + musclexposition(ch), electrodeyposition, 'k')
    plot( musclexposition(ch)*ones(1,6), linspace(1, 6, 6), '--k')
    %xlim([0 1])
end

sortedmuscles = {'DEL', 'BIC', 'FCR', 'TRI', 'EDC', 'ECR', 'FDS'}
title(pinNames{iF})
set(gca, 'ytick', [1:6], 'yticklabel', {'E6', 'E5', 'E4', 'E3', 'E2', 'E1'})
set(gca, 'xtick', [0.25:0.25:1.75], 'xticklabel',sortedmuscles)

legend(muscles ,'Location','NorthWest')
%% Triggered responses 
% All muscles are represented on the same scale one after the other, one subplot for each electrode
% The order of the muscles is 
% DEL BIC TRI EDC FCR ECR FDS ABP
h =figure('units','normalized','outerposition',[0 0 0.5 1]);
dims = numSubplots(length(pinNames));
for iF = 1 : length(fileNames) 
    subplot(dims(1), dims(2), iF)
    for iM = 1 : length(TDTdat{iF}.muscles)    

        
        plotEMGtrig_ontime(TDTdat{iF}, time(iM,:), TDTdat{iF}.muscles{iM}, muscleColors)
        ylim([MINACT MAXACT])
        title(pinNames{iF})
    end
   
%     saveas(h,...
%     fullfile('/Users/barrab/Documents/PhD/MkAnalysis/Figures/Brienne/20190603', ['TrigResponses_', fileNames{iF}]), ...
%     'png')
    
end

%% Plotting muscle triggered responses for each file, comparing electrodes for the same muscle

dims = numSubplots(length(pinNames));
maxactivation = 0 ; 
minactivation = 1000;
for iM = 1 : length(TDTdat{iF}.muscles)
   
    h =figure('units','normalized','outerposition',[0 0 0.5 1]);
    for iF = 1 : length(fileNames)  

        maxactivation = max([max(TDTdat{iF}.trigResp{iM}(:)), maxactivation]);
        minactivation = min([min(TDTdat{iF}.trigResp{iM}(:)), minactivation]);
    end
    for iF = 1 : length(fileNames)  

    
        subplot(dims(1), dims(2), iF)
        plotEMGtrig(TDTdat{iF}, TDTdat{iF}.muscles{iM}, muscleColors)
        ylim([minactivation maxactivation])
        title([pinNames{iF}, ' ',TDTdat{iF}.muscles{iM} ])
    end
%     saveas(h,...
%     fullfile('/Users/barrab/Documents/PhD/MkAnalysis/Figures/Brienne/20190603', ['TrigResponses_', fileNames{iF}]), ...
%     'png')
    
end


% %% Preparing data structs for Algorithm
% 
% [RClot] = computeRClot(TDTdat)
