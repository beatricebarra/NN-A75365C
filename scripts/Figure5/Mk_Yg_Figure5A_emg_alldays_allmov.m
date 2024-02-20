warning('off','all')
clear all
clc
close all
% Root = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/';
% currentpath = '/Users/Bea/Documents/PhD/Repos/unifr_cervicalEES/Code/scripts/';
Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES/Code/scripts/';


Animal = 'Ygritte';
expDates = {'20190809', '20190812','20190813','20190815','20190820',};
selectedFiles = {
    [4 5 6 11 14 16 17 ]; ... 
    [2 5 7 8 9 10 12 13 17 18]; ...
    [2 3   10  13 14]; ...
    [ 2 4 5 6  8]; ...
    [1 3 4 6 8 9 10 12 14 17]; ...
}% each cell is a date
stimMask = {
    [0 1 1  1 1  1 1]; ...
    [0 1 1 1 1 1 1 1 1 1 ]; ...
    [0 0  1 1 1 1 1]; ...
    [0 1 0 1 1]; ...
    [0 0 1 1 1 1 1 1 1 1]; ...
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
            TDTdat{iif}.iFile = iFile;
            
        end
        if any(strcmp(recsystems,'VICON'))
            VICONfile = [filenameVICON, num2str(iFile, '%02.f')];
            try 
                VICONdat{iif} = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '.mat' ]));
            catch
                VICONdat{iif} = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '_M.mat' ]));
            end
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
            TDTdat{iFile} = processEMGYgritte(TDTdat{iFile}, 'blankart', 'filter');%'blankart' 
        else
            TDTdat{iFile} = processEMGYgritte(TDTdat{iFile}, 'filter');%'blankart' 
        end

    end

    if any(strcmp(recsystems,'VICON'))
        %Computation of joint angles
        VICONdat{iFile}.jointAngles = computeJointAngles(VICONdat{iFile});

        % Retrieve analogs
        VICONdat{iFile} = getVICONAnalogIn(VICONdat{iFile}, analogs);


        % Syncronization
        % Number of square peaks in TDT trig in VICON
        VICONdat{iFile}.TDT_Trigger.times = alignTDT2VICON_robust(TDTdat{iFile}, VICONdat{iFile});
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
    iFile
    % Cleaning events
    [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(...
        TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}, ...
        {'startMov', 'grabMov', 'endMov' });
    % Extracting movement attempt indicator
    TDTdat{iFile} = extractMovementAttempt(TDTdat{iFile}, 0.5);
    % Transforming triggers in blackrock in burst triggers
    BKRdat{iFile} = extractBursts_BKR(BKRdat{iFile});
    % Clean force signals 
    VICONdat{iFile} = cleanForce(VICONdat{iFile});
    if TDTdat{iFile}.isstim == 1
        isstimtimed = checkStimTiming(BKRdat{iFile}); 
        TDTdat{iFile}.isstimtimed = isstimtimed;
        
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

            title(['Muscle ' , num2str(iM)] )
            bar(date_loc_plot(iD)+offset, mean(stim), 'r')
            errorbar(date_loc_plot(iD) +offset, mean(stim), 0, std(stim)) 
            randvect = 0.1*(rand(size(stim))-0.5); 
            plot(date_loc_plot(iD)*ones(size(randvect))+randvect+offset, stim, 'ob')
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


