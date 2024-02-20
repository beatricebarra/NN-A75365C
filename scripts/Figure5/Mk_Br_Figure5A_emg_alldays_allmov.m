warning('off','all')
clear all
clc
close all

Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES/Code/scripts/';
Animal = 'Brienne';


expDates = {  '20190608', '20190610', '20190614', '20190618', '20190624'};
selectedFiles = { 
    [1 2 3 4 5 6 7 8 9 10 11]; ... 
    [2 3  10 11 12 13 16 18 19]; ...
    [1 3 4 5 8 13];...
    [2, 3, 5, 6, 7, 8, 9];...
    [1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 20 ];...
}% each cell is a date
stimMask = {
    
    [0 0 0 1 1 1 1 1 0 0 1]; ...
    [0 0  1 1 1 1 1 1 1];...
    [0 0 1 1 1 0]; ... 
    [0 0 1 1 1 1 0];...
    [0 0 1 1 0 1 1 0 1 1 1 0 1 1 1 1 0 1]; ... 
}% each cell is a date, 1 means stim file, 0 is baseline


scspapercolors = [230 230 230; 125 125 125; 255 183 27];
recsystems = {'TDT', 'VICON', 'BLACKROCK'};
additionalpath = 'TDT';
load('/Users/barrab01/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')

set(0, 'DefaultFigureRenderer', 'painters');

iif = 0;
for iDate = 1 : length(selectedFiles)
    for iF = 1 : length(selectedFiles{iDate})
        iif = iif + 1;
        expDate = expDates{iDate};
        filenameTDT = [expDate,'_', Animal, '-'];
        if strcmp(expDate, '20190820')
            filenameTDT = [expDate,'_', Animal, 'HF-'];
        end
        if strcmp(expDate, '20190604')
           filenameTDT = [expDate,'_', Animal, '_-'];
        end
        filenameVICON = [expDate,'_',Animal, '_'];
        filenameBKR = [expDate,'_',Animal, '_'];
        dataPath = fullfile(Root, Animal,expDate, additionalpath);

        iFile = selectedFiles{iDate}(iF);
        disp(['Loading File n ', num2str(iFile) ])
        if any(strcmp(recsystems,'TDT')) || any(strcmp(recsystems,'TDTDyn'))
            TDTfile = [filenameTDT, num2str(iFile, '%01.f')];
            TDTdat{iif} = load(fullfile(dataPath, [TDTfile, '.mat'] ));
            TDTdat{iif}.fs = 12207.03;
            TDTdat{iif}.name = TDTfile;
            TDTdat{iif}.date = expDate;
            TDTdat{iif}.isstim = stimMask{iDate}(iF);
            
        end
        if any(strcmp(recsystems,'VICON'))
            VICONfile = [filenameVICON, num2str(iFile, '%02.f')];
            VICONdat{iif} = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '.mat' ]));
        end
        if any(strcmp(recsystems,'BLACKROCK'))
            BKRfile = [filenameBKR,num2str(iFile, '%03.f')];
            nsFile = fullfile(Root, Animal,expDate, 'BLACKROCK', [BKRfile, '.ns4'] ); 
            BKRdat{iif}.nsData = openNSx('read', nsFile,'p:double', 's:1');
           
        end
    end
    
end



%% EMG preprocessing stage:

% Extracting features from VICON data
% - joint angles
% - analog signals
% - events

% Extracting features from TDT data
% - bipolar EMGS
% - envelope

iif = 0;
analogs = {'BlackRock', 'TDT_Trigger', 'Kuka_Position'};
for iFile = 1 : length(TDTdat)
     
    
    disp(['Processing File n ', num2str(iFile) ])

    if any(strcmp(recsystems,'TDT'))
        % Preprocessing EMGs
        if TDTdat{iFile}.isstim == 1
            TDTdat{iFile} = processEMGBrienne(TDTdat{iFile}, 'blankart', 'filter');%'blankart' 
        else
            TDTdat{iFile} = processEMGBrienne(TDTdat{iFile}, 'filter');%'blankart' 
        end
        

    end

    if any(strcmp(recsystems,'VICON'))
        %Computation of joint angles
        VICONdat{iFile}.jointAngles = computeJointAngles(VICONdat{iFile});

        % Retrieve analogs
        VICONdat{iFile} = getVICONAnalogIn(VICONdat{iFile}, analogs);

        
        % Syncronization
         [VICONdat{iFile}] = alignTDT2VICON_allcases(TDTdat{iFile}, VICONdat{iFile}, Animal, TDTdat{iFile}.date);
        [peak, VICONdat{iFile}.BKR_trigger.samples] = (findpeaks(VICONdat{iFile}.analog.data(:,1), 'MinPeakHeight', 1, 'MinPeakProminence', 1));
        VICONdat{iFile}.BKR_trigger.times = VICONdat{iFile}.BKR_trigger.samples(1)/VICONdat{iFile}.analog.framerate;

        % Retrieve events
        VICONdat{iFile} = getVICONEvents(VICONdat{iFile});
        
        
    end
    
    if any(strcmp(recsystems,'BLACKROCK'))
        BKRdat{iFile}.TDT_Trigger.times = VICONdat{iFile}.TDT_Trigger.times - VICONdat{iFile}.BKR_trigger.times;
        BKRdat{iFile}.TDT_Trigger.samples = floor(BKRdat{iFile}.TDT_Trigger.times * BKRdat{iFile}.nsData.MetaTags.SamplingFreq);
    end
end

%% Additional processing (more advanced and optional)

for iFile = 1 : length(TDTdat)
    deltaBKRTDT = VICONdat{iFile}.TDT_Trigger.times - VICONdat{iFile}.BKR_trigger.times; 
    % Cleaning events
    [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(...
        TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}, ...
        {'startMov', 'grabMov', 'endMov'});
    % Extracting movement attempt indicator
    %TDTdat{iFile} = extractMovementAttempt(TDTdat{iFile}, 0.5, [1 2], 0.15);
    % Transforming triggers in blackrock in burst triggers
    BKRdat{iFile} = extractBursts_BKR(BKRdat{iFile});
    % Clean BKR bursts
    BKRdat{iFile} = cleanBursts_BKR(BKRdat{iFile}, Animal, TDTdat{iFile}.date);
    % Clean force signals 
    VICONdat{iFile} = cleanForceMkBr(VICONdat{iFile});
    if TDTdat{iFile}.isstim == 1
        isstimtimed = checkStimTiming(BKRdat{iFile}); 
        TDTdat{iFile}.isstimtimed = isstimtimed;
        BKRdat{iFile}.reachstim_sniptime = extractReachStim(BKRdat{iFile}); 
        BKRdat{iFile}.pullstim_sniptime = extractPullStim(BKRdat{iFile}); 
        for iE = 1 : length(BKRdat{iFile}.reachstim_sniptime)
            TDTdat{iFile}.reachstim_sniptime{iE} = BKRdat{iFile}.reachstim_sniptime{iE} - deltaBKRTDT;
            TDTdat{iFile}.pullstim_sniptime{iE} = BKRdat{iFile}.pullstim_sniptime{iE} - deltaBKRTDT; 
        end
    end

end


%% STIM effect on the whole movement 
% Available Pins : 1 2 and 5

dt = 1/TDTdat{1}.fs;
NOSTIM = cell(length(expDates), 1);
for i = 1 : length(NOSTIM)
    NOSTIM{i} = [];
end
STIM = cell(length(expDates), 1);
for i = 1 : size(STIM, 1)
     STIM{i} = [];
end
    
for iFile = 1 : length(TDTdat)% [ 1:22 , 30:33, 34:length(TDTdat) ]
    

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
               
               
                if TDTdat{iFile}.isstimtimed.reach(iM) ==1 || TDTdat{iFile}.isstimtimed.pull(iM) ==1
                    for iC = 1 : size(TDTdat{iFile}.EMGb,2)
                        temp(iC) = sum(dt.*abs(TDTdat{iFile}.EMGb(startM:endM, iC)).^2)...
                                    ./length(TDTdat{iFile}.EMGb(startM:endM, iC));
                    end
                    STIM{date_idx} = [STIM{date_idx}; temp];
                end
                            
            end
        end
    end
    
end

% Plot

figure

date_loc_plot = [2 4 6 8 10]; 

offset = 0.5; 
npoints= cell(8,  size(NOSTIM,1));
for iM =1:8
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
            npoints{iM, iD} = [npoints{iM, iD}; length(baseline)];

            title(['Muscle ' , num2str(iM), ', Pin ', num2str(iE)])
            bar(date_loc_plot(iD)+offset, mean(stim), 'r')
            errorbar(date_loc_plot(iD) +offset, mean(stim), 0, std(stim)) 
            randvect = 0.1*(rand(size(stim))-0.5); 
            plot(date_loc_plot(iD)*ones(size(randvect))+randvect+offset, stim, 'ob')
            npoints{iM, iD} = [npoints{iM, iD}; length(stim)];
            p = ranksum(baseline, stim); 

            disp(['Muscle', TDTdat{iFile}.muscles{iM}, ' Date ', expDates{iD}, ' pvalue=' , num2str(p)])
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


