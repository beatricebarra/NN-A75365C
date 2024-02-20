%%
% @Author: barrab
% @Date:   2019-06-13 12:14:23
% @Last Modified by:   barrab
% @Last Modified time: 2019-06-13 12:14:23
% run with 
% 
% options.format = 'pdf';
% options.showCode = false;
% publish('Brienne',options);
%% Recruitment curves have been performed on : 
% 
% * 03.06.2019
% * 12.06.2019
% * 19.06.2019
% * 02.07.2019

% Relative parameters are 

Animal = 'Brienne';
pittPC =1;

params_12062019;
params_19062019;
params_02072019;
params_03062019;


params_19062019;

nResp = 4;
nChan = 8;
fS = 12207.03;
lengthResp = floor(0.05*fS);
singleLead = 0;
OneTrueCM;
emgmapping_brienne;
% Loading files and preprocessing

for i = 1 : length(fileNames)
    
    % Loading
    TDTdat{i} = load(fullfile(Root, [fileNames{i}, '.mat']));
    % Saving additional infos
    TDTdat{i}.currents = currents{i};
    TDTdat{i}.muscles = muscles;
    
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
    
    % Number of explored currents during recruitment curves
    nCurrents = size(TDTdat{i}.EMGtb, 1)/(lengthResp * nResp);
    % Computing recruitment curves based on peak to peak measure
    for ch = 1 : nChan
        for iCurr = 1 : nCurrents
            allResps = reshape(TDTdat{i}.EMGtb((iCurr-1)*lengthResp*nResp+1: (iCurr)*lengthResp*nResp, ch), [lengthResp, nResp]);
            allResps(1:42,:) = 0;
            avgResp = median(allResps, 2);
            TDTdat{i}.trigResp{ch}(:,:,iCurr) = allResps;
            TDTdat{i}.p2p(ch, iCurr) = max(avgResp) - min(avgResp);
        end
    end
    
    
end

% for iF = 1 : length(fileNames)
%     TDTdat{iF}.currents = currents{iF};
%     TDTdat{iF}.electrode = pinNames{iF};
% end
% time = linspace(0, 0.05, 610);
% for i = 1 : length(fileNames)
%     disp(fileNames{i})
%     for ch = 7
%         nCurrents = size(TDTdat{i}.EMGtb, 1)/(lengthResp * nResp);
%         for iCurr = 1 : nCurrents
%         plot(time, TDTdat{i}.trigResp{ch}(:,:,iCurr))
%         pause
%         clf
%         end
%     end
% end

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
%save('/Users/Bea/Documents/PhD/Repos/unifr_cervicalEES/Structs/RCBrienne', 'RCdat')

%% Computing selectivity index
clear allRC
for i = 1 : length(fileNames)
    allRC{i}.table = TDTdat{i}.RC;
end
SI_total = average_SI(allRC, length(pinNames), length(muscles));
figure
deltatheta = 360/size(SI_total, 1);
theta = (linspace(deltatheta/2, 360-deltatheta/2, size(SI_total, 1))./360*2*pi);
for iP = 1 : size(SI_total, 2)
    subplot(2,4,iP)
    ps = polarplot([theta theta(1)],[SI_total(:,iP) ;SI_total(1,iP)], 'k')
    hold on
    for iM = 1 : length(muscles)
        ps = polarscatter(theta(iM),SI_total(iM,iP))
        ps.SizeData = 40;
        ps.MarkerEdgeColor = muscleColors.(muscles{iM});
        ps.MarkerFaceColor = muscleColors.(muscles{iM});
    end
    rlim([0 0.4])
    
    thetaticks(theta*360/(2*pi))
    thetaticklabels(muscles)
    %set(gca, 'thetatick', theta, 'thetaticklabels', muscles)
end

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



%% Plotting activation
rostrocaudalindex = [8 7 6 3 4 5 2 1];

h =figure('units','normalized','outerposition',[0 0 0.5 1]);
dims = numSubplots(length(fileNames) );

for iF = 1 : length(fileNames) 
    subplot(dims(1), dims(2), iF)
    hold on
    for ch = 1 : nChan
        activationindex(ch) = mean(TDTdat{iF}.RC(ch,:));
        plot(activationindex(ch), rostrocaudalindex(ch),  'o', 'Linewidth', 2,...
            'MarkerFaceColor', muscleColors.(muscles{ch}), ...
            'MarkerEdgeColor', muscleColors.(muscles{ch}), ...
            'Markersize', 10)
    end
    [sortrcindex, isort] = sort(rostrocaudalindex);
    plot(activationindex(isort), sortrcindex, 'k')
    xlim([0 1])
end


title(pinNames{iF})
legend(muscles ,'Location','NorthWest')


%% Plotting activation index as the motor pools

rostrocaudalindex = [8 7 6 3 4 5 2 1];
musclexposition = [1 2 4 5 3 6 7 8]/4;
    %]/2;
electrodeyposition = [6 5 4 3 2 1];
h =figure('units','normalized','outerposition',[0 0 0.5 1]);
hold on

for ch = 1 : nChan-1
    for iF = 1 : length(fileNames) 
        activationindex(iF) = mean(TDTdat{iF}.RC(ch,:));
        plot(activationindex(iF) + musclexposition(ch), electrodeyposition(iF),  'o', 'Linewidth', 2,...
            'MarkerFaceColor', muscleColors.(muscles{ch}), ...
            'MarkerEdgeColor', muscleColors.(muscles{ch}), ...
            'Markersize', 10)
    end
    
    plot(activationindex + musclexposition(ch), electrodeyposition, 'k')
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




%%


% %% Plotting muscle triggered responses for each file, each muscle each current
% 
% dims = numSubplots(8);
% for iF = 1 : length(fileNames)  
%     
%     h =figure('units','normalized','outerposition',[0 0 1 1])
%     for iM = 1 : length(TDTdat{iF}.muscles)
%         subplot(dims(1), dims(2), iM)
%         plotEMGtrig(TDTdat{iF}, TDTdat{iF}.muscles{iM}, muscleColors)
%         title(TDTdat{iF}.muscles{iM})
%     end
% %     saveas(h,...
% %     fullfile('/Users/barrab/Documents/PhD/MkAnalysis/Figures/Brienne/20190603', ['TrigResponses_', fileNames{iF}]), ...
% %     'png')
%     
% end



% Preparing data structs for Algorithm
%[RClot] = computeRClot(TDTdat)
