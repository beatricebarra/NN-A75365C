warning('off','all')
clear all
clc
close all

Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES/Code/scripts/NatNeuro_finalsubmission';
load('/Users/barrab01/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')

Animal = 'Ygritte';

expDates = { '20190813', '20190815' };%

selectedFiles = {
    
   %[ 4 5 10 11 12 14 15 16 ];...% 09-08
   %[ 2  5 13  17 18]; ...% 12-08
   [ 11 12 13 14 15 16]; ...%13-08
   [3 4  6 7 8] ; ...% 15-08
   %[1 3 4 6 8 9 10 11 12 14 17]; ...%20 -08
    }
stimMask = {
    
    %[0 0 1 1 1 1 1 1]; ...
    %[ 0 1 1 1 1]
    [  1 1 1 1 1 1]; ...%0
    [ 1 1  1 1 1 ]; 
    %[0 0 1 1 1 1 1 1 1 1 1]; ...%0  

}% each cell is a date, 1 means stim file, 0 is baseline

scspapercolors = [230 230 230; 125 125 125; 255 183 27];
recsystems = {'TDT', 'VICON', 'BLACKROCK'};
additionalpath = 'TDT';

set(0, 'DefaultFigureRenderer', 'painters');

iif = 0;
for iDate = 1 : length(selectedFiles)
    [selected_matfiles] = getNames_reachEventFiles(Animal, expDates{iDate});
    eventfiles_numbers = cellfun(@(x) str2num(x(18:19)),selected_matfiles );
    for iF = 1 : length(selectedFiles{iDate})
        iif = iif + 1;
        expDate = expDates{iDate};
        filenameTDT = [expDate,'_', Animal, '-'];
        if strcmp(expDate, '20190820')
            filenameTDT = [expDate,'_', Animal, 'HF-'];
        end
        if strcmp(expDate, '20190529')|| strcmp(expDate, '20190nsecs4')
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
            idx_eventfile = find(eventfiles_numbers == iFile);
            VICONfile = [filenameVICON, num2str(iFile, '%02.f')];
            temp = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '_kin.mat' ]));
            VICONdat{iif}.kinematic = temp.kinematic;
            if stimMask{iDate}(iF) == 2
                %temp = load(fullfile(Root, Animal,expDate, 'VICON', [eventfilename]));
                temp = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '.mat']));
                VICONdat{iif}.event = temp.event;
                temp = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '.mat']));
                VICONdat{iif}.analog = temp.analog;
            else
                temp = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '.mat']));
                VICONdat{iif}.event = temp.event;
                VICONdat{iif}.analog = temp.analog;
            end
           
            
        end
        if any(strcmp(recsystems,'BLACKROCK'))
            BKRfile = [filenameBKR,num2str(iFile, '%03.f')];
            nsFile = fullfile(Root, Animal,expDate, 'BLACKROCK', [BKRfile, '.ns4'] );
            %nevfile = fullfile(Root, Animal,expDate, 'BLACKROCK', [BKRfile, '.nev'] );
            BKRdat{iif}.nsData = openNSx('read', nsFile,'p:double', 's:1');
            %BKRdat{iif}.nevData = openNEV(nevfile,'nomat','nosave');
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
        %figh = plotEMG(TDTdat{iFile}, 'envonraw')
        %figh = plotEMG(TDTdat{iFile}, 'envonrawn')

    end

    if any(strcmp(recsystems,'VICON'))
        %Computation of joint angles
        VICONdat{iFile}.jointAngles = computeDLCJointAngles(VICONdat{iFile});

        % Retrieve analogs
        VICONdat{iFile} = getVICONAnalogIn(VICONdat{iFile}, analogs);

        %Retrieve KUKA position
        %VICONdat{iFile} = getKUKAposition(VICONdat{iFile});

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
    % Cleaning events
    [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(...
        TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}, ...
        {'startMov', 'grabMov', 'endMov' });
    % Extracting movement attempt indicator
    % TDTdat{iFile} = extractMovementAttempt(TDTdat{iFile}, 0.5);
    % Transforming triggers in blackrock in burst triggers
    BKRdat{iFile} = extractBursts_BKR(BKRdat{iFile});
    % Clean force signals 
    %VICONdat{iFile} = cleanForce(VICONdat{iFile});
    VICONdat{iFile} = retrieveBadTrials(VICONdat{iFile}); 
end


%% Check one trial at a time
% BAS = 1 4 7
% STIM  2 3 5 6 8 9 10
for iFile = 1
    figure
    completetrialPlot(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}, {'EMG', 'stim'})
    
end

%% determine which trials are early stim

for iFile = 1 : length(VICONdat)
    startMovs = VICONdat{iFile}.events_bad.startMov; 
    endMovs = VICONdat{iFile}.events_bad.endMov; 
    early_label = zeros(size(startMovs)); 
    for iM = 1 : length(startMovs)
        idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        stimstarts = BKRdat{iFile}.burstsStart{6} + VICONdat{iFile}.BKR_trigger.times; 
        idxburst = intersect(find(stimstarts > startMovs(iM)), find(stimstarts < endMovs(iM))); 
        if not(isempty(idxburst))
            timestart = stimstarts(idxburst(1)); 
            idx_timestart = floor(timestart*VICONdat{iFile}.kinematic.framerate); 
            if  -VICONdat{iFile}.kinematic.x(idx_timestart,end)< -80
                early_label(iM) = 1;
                early_stim_time(iM) = timestart; 
            end
        end
           
    end
    VICONdat{iFile}.events_bad.early_label = early_label;
    VICONdat{iFile}.events_bad.early_stim_time = early_stim_time;
end


for iFile = 1 : length(VICONdat)
    startMovs = VICONdat{iFile}.events_tokeep.startMov; 
    endMovs = VICONdat{iFile}.events_tokeep.endMov; 
    timestart = zeros(size(startMovs)); 
    for iM = 1 : length(startMovs)
        idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        stimstarts = BKRdat{iFile}.burstsStart{6} + VICONdat{iFile}.BKR_trigger.times; 
        idxburst = intersect(find(stimstarts > startMovs(iM)), find(stimstarts < endMovs(iM))); 
        if not(isempty(idxburst))
            timestart(iM) = stimstarts(idxburst(1)); 
        end
           
    end
    
    VICONdat{iFile}.events_tokeep.stim_time = timestart;
end

%% Plot average trajectories for bad and good
figure
hold on

allx = [];
ally = [];
allz = [];

for iFile = 1 : length(VICONdat)
    startMovs = VICONdat{iFile}.events_bad.startMov; 
    endMovs = VICONdat{iFile}.events_bad.endMov; 
    for iM = 1 : length(startMovs)
        if VICONdat{iFile}.events_bad.early_label(iM) == 1
            idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
            idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 

            x = VICONdat{iFile}.kinematic.x(idxStart:idxStop,end); 
            y = VICONdat{iFile}.kinematic.y(idxStart:idxStop,end); 
            z = VICONdat{iFile}.kinematic.z(idxStart:idxStop,end); 
            t_old = linspace(1, length(x), length(x)); 
            old_idx_stim = floor(VICONdat{iFile}.events_bad.early_stim_time(iM)*VICONdat{iFile}.kinematic.framerate) - idxStart;
            t_new = linspace(1, length(x), 1000); 
            new_idx_stim = floor(old_idx_stim/ length(t_old)*1000); 
            x = interp1(t_old, x, t_new); 
            y = interp1(t_old, y, t_new); 
            z = interp1(t_old, z, t_new);


            allx = [allx; -x]; 
            ally = [ally; y]; 
            allz = [allz; z]; 
            plot(new_idx_stim, -x(new_idx_stim), 'or')
        end
        %plot(x, 'r')
    end
    
end

stdshade(allx, 0.2, 'r')
%
allx = [];
ally = [];
allz = [];
for iFile = 1 : length(VICONdat)
    startMovs = VICONdat{iFile}.events_tokeep.startMov
    endMovs = VICONdat{iFile}.events_tokeep.endMov; 
    for iM = 1 : length(startMovs)
        idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        
        x = VICONdat{iFile}.kinematic.x(idxStart:idxStop,end); 
        y = VICONdat{iFile}.kinematic.y(idxStart:idxStop,end); 
        z = VICONdat{iFile}.kinematic.z(idxStart:idxStop,end); 
        t_old = linspace(1, length(x), length(x)); 
        t_new = linspace(1, length(x), 1000); 
        x = interp1(t_old, x, t_new); 
        y = interp1(t_old, y, t_new); 
        z = interp1(t_old, z, t_new);
        
       
        allx = [allx; -x]; 
        ally = [ally; y]; 
        allz = [allz; z]; 
        %plot(x, 'k')
    end
end

stdshade(allx, 0.2, 'k')

%% Find normalization value
%% Plot trajectory good trials with stimulation marker

max_reached = [];
min_reached = [];
for iFile = 2 : 7%length(VICONdat)
    startMovs = VICONdat{iFile}.events_bad.startMov; 
    endMovs = VICONdat{iFile}.events_bad.startMov + 1; 
    for iM = 1 : length(startMovs)
        if VICONdat{iFile}.events_bad.early_label(iM) == 1
            idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
            idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 

            x = VICONdat{iFile}.kinematic.x(idxStart:idxStop,end); 
            x = -x + x(1); 
            max_reached = [max_reached, max(x)]; 
           
        end
    end
end

for iFile = 1 : length(VICONdat)
    startMovs = VICONdat{iFile}.events_tokeep.startMov
    endMovs = VICONdat{iFile}.events_tokeep.startMov + 1; 
    for iM = 1 : length(startMovs)
        idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        
        x = VICONdat{iFile}.kinematic.x(idxStart:idxStop,end); 
        x = -x + x(1); 
        max_reached = [max_reached, max(x)]; 
        
       
    end
end


max_norm = max(max_reached); 
%% Plot trajectory good trials with stimulation marker
figure
hold on

allx = [];
ally = [];
allz = [];

for iFile = 2 : 7%length(VICONdat)
    startMovs = VICONdat{iFile}.events_bad.startMov; 
    endMovs = VICONdat{iFile}.events_bad.startMov + 1; 
    for iM = 1 : length(startMovs)
        if VICONdat{iFile}.events_bad.early_label(iM) == 1
            idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
            idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 

            x = VICONdat{iFile}.kinematic.x(idxStart:idxStop,end); 
            x = -x + x(1);
            x = x / max_norm; 
            y = VICONdat{iFile}.kinematic.y(idxStart:idxStop,end); 
            z = VICONdat{iFile}.kinematic.z(idxStart:idxStop,end); 
            t_old = linspace(1, length(x), length(x)); 
            old_idx_stim = floor(VICONdat{iFile}.events_bad.early_stim_time(iM)*VICONdat{iFile}.kinematic.framerate) - idxStart;
%             t_new = linspace(1, length(x), 1000); 
%             new_idx_stim = floor(old_idx_stim/ length(t_old)*1000); 
%             x = interp1(t_old, x, t_new); 
%             y = interp1(t_old, y, t_new); 
%             z = interp1(t_old, z, t_new);

            if old_idx_stim < (idxStop - idxStart)
                allx = [allx; x']; 
                ally = [ally; y]; 
                allz = [allz; z]; 
                plot(old_idx_stim, x(old_idx_stim), 'or')
                plot(t_old, x, 'r')
                iFile
                iM
                pause
            end
        end
        %plot(x, 'r')
    end
    
end
%stdshade(allx, 0.2, 'r')
hold on
allx = [];
ally = [];
allz = [];
for iFile = 1 : length(VICONdat)
    startMovs = VICONdat{iFile}.events_tokeep.startMov
    endMovs = VICONdat{iFile}.events_tokeep.startMov + 1; 
    for iM = 1 : length(startMovs)
        idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        
        x = VICONdat{iFile}.kinematic.x(idxStart:idxStop,end); 
        x = -x + x(1); 
        x = x / max_norm; 
        y = VICONdat{iFile}.kinematic.y(idxStart:idxStop,end); 
        z = VICONdat{iFile}.kinematic.z(idxStart:idxStop,end); 
        t_old = linspace(1, length(x), length(x)); 
        old_idx_stim = floor(VICONdat{iFile}.events_tokeep.stim_time(iM)*VICONdat{iFile}.kinematic.framerate) - idxStart;
        if old_idx_stim> (idxStop-idxStart)
            %plot(t_old, -x, 'k')
        elseif old_idx_stim> 0
            plot(old_idx_stim, x(old_idx_stim), 'ok')
            %plot(t_old, -x, 'k')
        end
%         t_new = linspace(1, length(x), 1000); 
%         x = interp1(t_old, x, t_new); 
%         y = interp1(t_old, y, t_new); 
%         z = interp1(t_old, z, t_new);
        
       
        allx = [allx; x']; 
        ally = [ally; y']; 
        allz = [allz; z']; 
       
    end
end
hold on
stdshade(allx, 0.2, 'k')



%% Max x reached

figure
hold on

max_reachx_early = [];
for iFile = 2 : 7%length(VICONdat)
    startMovs = VICONdat{iFile}.events_bad.startMov; 
    endMovs = VICONdat{iFile}.events_bad.startMov + 1; 
    for iM = 1 : length(startMovs)
        if VICONdat{iFile}.events_bad.early_label(iM) == 1
            idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
            idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 

            x = VICONdat{iFile}.kinematic.x(idxStart:idxStop,end)- VICONdat{iFile}.kinematic.x(idxStart, end); 
            x = -x + x(1); 
            x = x / max_norm; 
            old_idx_stim = floor(VICONdat{iFile}.events_bad.early_stim_time(iM)*VICONdat{iFile}.kinematic.framerate) - idxStart;

            if old_idx_stim < (idxStop - idxStart)
                iFile, iM
                pause
                max_reachx_early = [max_reachx_early; max(x)] ;    
            end
        end
    end
    
end
max_reachx_good = []
for iFile = 1 : length(VICONdat)
    startMovs = VICONdat{iFile}.events_tokeep.startMov
    endMovs = VICONdat{iFile}.events_tokeep.startMov + 1; 
    for iM = 1 : length(startMovs)
        idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
        
        x = VICONdat{iFile}.kinematic.x(idxStart:idxStop,end)- VICONdat{iFile}.kinematic.x(idxStart, end); 
        x = -x + x(1); 
        x = x / max_norm; 
        old_idx_stim = floor(VICONdat{iFile}.events_tokeep.stim_time(iM)*VICONdat{iFile}.kinematic.framerate) - idxStart;
        if old_idx_stim> (idxStop-idxStart)
            max_reachx_good = [max_reachx_good; max(x)] ;
            iFile, iM
        elseif old_idx_stim> 0
            max_reachx_good = [max_reachx_good; max(x)] ;  
            iFile, iM
        end       
    end
end
figure
hold on
plot(ones(size(max_reachx_early)),max_reachx_early, 'o' )
bar(1, mean(max_reachx_early)); 
errorbar(1, mean(max_reachx_early), 0, std(max_reachx_early)); 
plot(2*ones(size(max_reachx_good)),max_reachx_good, 'o' )
bar(2, mean(max_reachx_good)); 
errorbar(2, mean(max_reachx_good), 0, std(max_reachx_good)); 
p = ranksum(max_reachx_early,max_reachx_good ); 
title(num2str(p))
%%
% plot one example stick diagram
iFile = 4;
startMovs = VICONdat{iFile}.events_tokeep.startMov; 
endMovs = VICONdat{iFile}.events_tokeep.startMov + 1; 
iM = 1; 
idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
shoulder.x = VICONdat{iFile}.kinematic.x(:,1);  
shoulder.y= VICONdat{iFile}.kinematic.z(:,1); 
elbow.x = VICONdat{iFile}.kinematic.x(:,2);  
elbow.y= VICONdat{iFile}.kinematic.z(:,2); 
wrist.x = VICONdat{iFile}.kinematic.x(:,3);  
wrist.y= VICONdat{iFile}.kinematic.z(:,3); 
figure
stickdiagram_ygritte_behav(shoulder, elbow, wrist, [idxStart, idxStop], 0, 1, [0 0 0])


%%
hold on
iFile = 7;
startMovs = VICONdat{iFile}.events_bad.startMov; 
endMovs = VICONdat{iFile}.events_bad.startMov + 1; 
iM =1; 
idxStart = floor(startMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
idxStop = floor(endMovs(iM)*VICONdat{iFile}.kinematic.framerate); 
shoulder.x = VICONdat{iFile}.kinematic.x(:,1);  
shoulder.y= VICONdat{iFile}.kinematic.z(:,1); 
elbow.x = VICONdat{iFile}.kinematic.x(:,2);  
elbow.y= VICONdat{iFile}.kinematic.z(:,2); 
wrist.x = VICONdat{iFile}.kinematic.x(:,3);  
wrist.y= VICONdat{iFile}.kinematic.z(:,3); 

stickdiagram_ygritte_behav(shoulder, elbow, wrist, [idxStart, idxStop], 0, 1, [1 0 0])