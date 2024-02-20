clc
clear all 
close all
warning('off','all')

%% 

% Root = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/';
% currentpath = '/Users/Bea/Documents/PhD/Repos/unifr_cervicalEES/Code/scripts/';
Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
%currentpath = 'C:\GIT\unifr_cervicalEES\';

Animal = 'Sansa';


expDates = {  '20180924'};
filemap_20190924;
selectedFiles = {
    [ 1: 10, 12: 19 ]; ... [1 7 15 16 19 ]
}% each cell is a date
stimMask = {
    [0 1 1 1 1 1 0 1 1  1 1 1 1 1 1 1 1 1]; ... 
}% each cell is a date, 1 means stim file, 0 is baseline


scspapercolors = [230 230 230; 125 125 125; 255 183 27];
recsystems = {'TDT', 'VICON', 'BLACKROCK'};
additionalpath = 'TDT';
set(0, 'DefaultFigureRenderer', 'painters');

iif = 0;
for iDate = 1 : length(selectedFiles)
    for iF = 1 : length(selectedFiles{iDate})
        iif = iif + 1;
        expDate = expDates{iDate};
        filenameTDT = [];
        if strcmp(expDate, '20190820')
            filenameTDT = [expDate,'_', Animal, 'HF-'];
        end
        filenameVICON = [expDate,'_Sansa_Brain-'];
        filenameBKR = [expDate,'_Sansa_Brain'];% Example 20180413_KRG_TonicStim002.ns5
        dataPath = fullfile(Root, Animal,expDate, additionalpath);
        iFile = selectedFiles{iDate}(iF);
        disp(['Loading File n ', num2str(iFile) ])
        if any(strcmp(recsystems,'TDT')) || any(strcmp(recsystems,'TDTDyn'))
            TDTfile = [fileNames{iFile}];
            TDTdat{iif} = load(fullfile(dataPath, [TDTfile] ));
            TDTdat{iif}.fs = 12207.03;
            TDTdat{iif}.name = TDTfile;
            TDTdat{iif}.date = expDate;
            TDTdat{iif}.isstim = stimMask{iDate}(iF);
            
        end
        if any(strcmp(recsystems,'VICON'))
            VICONfile = [filenameVICON, num2str(iFile, '%02.f')];
            VICONdat{iif} = load(fullfile(Root, Animal, expDate,  'VICON', [VICONfile, '.mat' ]));
        end
        if any(strcmp(recsystems,'BLACKROCK'))
            BKRfile = [filenameBKR,num2str(iFile, '%03.f')];
            nsFile = fullfile(Root, Animal,expDate, 'BLACKROCK', [BKRfile, '.ns5'] );
            nevfile = fullfile(Root, Animal,expDate, 'BLACKROCK', [BKRfile, '.nev'] );
            %BKRdat{iif}.nsData = openNSx('read', nsFile,'p:double', 's:1');
            BKRdat{iif}.nevData = openNEV(nevfile,'nomat','nosave');
        end
    end
    
end


%% Additional processing



for iFile = [1 7 14 15 18]
    
    % Extracting comments from Blackrock --> in this case this is useful 
    % for the TDT trigger
    correction = 0;
    BKRdat{iFile} = extractComments_fromBKR(BKRdat{iFile} , correction);
    
    % VICON
    VICONdat{iFile} = getVICONEvents(VICONdat{iFile});
    
    LPbkrtrig = (VICONdat{iFile}.analog.data(:,1) >2);
    [peak, VICONdat{iFile}.BKR_Trigger.samples] = findpeaks(double(LPbkrtrig), 'MinPeakProminence', 0.8);
    VICONdat{iFile}.BKR_Trigger.times = VICONdat{iFile}.BKR_Trigger.samples/1000;
    VICONdat{iFile}.TDT_Trigger.times = VICONdat{iFile}.BKR_Trigger.times + BKRdat{iFile}.TDT_Trigger.times;
    
    % TDT
    TDTdat{iFile} = processEMGSansa(TDTdat{iFile}, 'filter', 'bipolar');
    
    % All together
    %[TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = copyEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile});
    [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = keepFullMov(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile});

end
%%
% Finding mean duration of movement to break down baseline trials
durations = [];
for iFile =[1 7 14 15 18]
    for iM = 1 : length(TDTdat{iFile}.events_tokeep.endMov)
        durations = [durations; ...
            TDTdat{iFile}.events_tokeep.endMov(iM) - TDTdat{iFile}.events_tokeep.startMov(iM)];
    end
end


%% Energy computation and plot
dt = 1/TDTdat{1}.fs;
NOSTIM = cell(length(expDates), 1);
for i = 1 : length(NOSTIM)
    NOSTIM{i} = [];
end
STIM = cell(length(expDates), 1);
for i = 1 : size(STIM, 1)
     STIM{i} = [];
end
    
for iFile = [1 7 14 15 18]% [ 1:22 , 30:33, 34:length(TDTdat) ]
    

    this_date = TDTdat{iFile}.date; 
    date_idx = find(cellfun(@(x) strcmp(x,this_date), expDates,'UniformOutput', true)); 
    for iM = 1 : length(TDTdat{iFile}.events_tokeep.startMov)
        startM = floor(TDTdat{iFile}.events_tokeep.startMov(iM)*TDTdat{iFile}.fs);
        endM = floor(TDTdat{iFile}.events_tokeep.endMov(iM)*TDTdat{iFile}.fs);
        
        if startM > 0 && endM <  size(TDTdat{iFile}.EMGb, 1)
            
            if TDTdat{iFile}.isstim == 0 % non stim files
                
                for iC = 1 : size(TDTdat{iFile}.EMGb,2)
                    temp(iC) = sum(dt.*abs(TDTdat{iFile}.EMGb(startM:endM, iC)).^2)...
                                ./length(TDTdat{iFile}.EMGb(startM:endM, iC));
                end
                NOSTIM{date_idx} = [NOSTIM{date_idx}; temp]; 

            else % stim files
               
               
                
                for iC = 1 : size(TDTdat{iFile}.EMGb,2)
                    temp(iC) = sum(dt.*abs(TDTdat{iFile}.EMGb(startM:endM, iC)).^2)...
                                ./length(TDTdat{iFile}.EMGb(startM:endM, iC));
                end
                STIM{date_idx} = [STIM{date_idx}; temp];
                
                            
            end
        end
    end
    
end

% Plot

figure

date_loc_plot = [2 4 6 8 10]; 

offset = 0.5; 
for iM =1:3
    for iD = 1 : size(NOSTIM,1)
        
        subplot(2, 4, iM)
        hold on

        if not(isempty(STIM{iD}))
            baseline = NOSTIM{iD}(:,iM);
            stim = STIM{iD}(:,iM); 
            bar(date_loc_plot(iD)-offset, mean(baseline), 'k')
            errorbar(date_loc_plot(iD)-offset, mean(baseline), 0, std(baseline))
            randvect = 0.1*(rand(size(baseline))-0.5); 
            plot(date_loc_plot(iD)*ones(size(randvect))+randvect- offset, baseline, 'og')

            title(['Muscle ' , num2str(iM)] )
            bar(date_loc_plot(iD)+offset, mean(stim), 'r')
            errorbar(date_loc_plot(iD) +offset, mean(stim), 0, std(stim)) 
            randvect = 0.1*(rand(size(stim))-0.5); 
            plot(date_loc_plot(iD)*ones(size(randvect))+randvect+offset, stim, 'ob')
            p = ranksum(baseline, stim)
            if p <0.05
                text(date_loc_plot(iD), mean([baseline; stim]), '*')
            end
        else 
            baseline = [];
            stim = [];
        end

        %ylim([0 5*10^-12])
        
    end
end
