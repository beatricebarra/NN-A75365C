warning('off','all')
clear all
clc
close all

Root = fullfile('/Volumes/SanDisk_Bea/SCS_Monkeys/');
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES/';
Animal = 'Brienne';
%expDates = {'20190608', '20190611',  '20190614', '20190624'};
expDates = { '20190529'};
selectedFiles = {
    [5:8 14:17]; ...This brings problems
}% each cell is a date
stimMask = {
    [2 2 2 2 2 2 2 2 ]; ...
}% each cell is a date, 1 means stim file, 0 is baseline
% expDates = { '20190529'};
% selectedFiles = {
%     [5:6]; ...This brings problems
% }% each cell is a date
% stimMask = {
%     [2 2 ]; ...
% }% each cell is a date, 1 means stim file, 0 is baseline


scspapercolors = [230 230 230; 125 125 125; 255 183 27];
recsystems = {'TDT', 'VICON', 'BLACKROCK'};
additionalpath = 'TDT';
brain_Fs = 30000;
artremovalparams.num_chan = 30;
artremovalparams.rejection_window = 0.0005;
artremovalparams.remove_flag =true;
%load('/Users/Bea/Documents/PhD/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')
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
        if strcmp(expDate, '20190529')
            filenameTDT = [expDate,'_', Animal, '_-'];
        end
        filenameVICON = [expDate,'_',Animal, '_'];
        filenameBKR = [expDate,'_',Animal, '_'];
        dataPath = fullfile(Root, Animal,expDate, additionalpath);

        iFile = selectedFiles{iDate}(iF);
        
        disp(['Loading File n ', num2str(iFile) ])
        if any(strcmp(recsystems,'TDT')) || any(strcmp(recsystems,'TDTDyn'))
            TDTfile = [filenameTDT, num2str(iFile-1, '%01.f')];
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
            nevfile = fullfile(Root, Animal,expDate, 'BLACKROCK', [BKRfile, '.nev'] );
            BKRdat{iif}.nsData = openNSx('read', nsFile,'p:double', 's:1');
            BKRdat{iif}.nevData = openNEV(nevfile,'nomat','nosave');
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
        %figh = plotEMG(TDTdat{iFile}, 'envonraw')
        %figh = plotEMG(TDTdat{iFile}, 'envonrawn')

    end

    if any(strcmp(recsystems,'VICON'))
        %Computation of joint angles
         VICONdat{iFile}.jointAngles = computeJointAngles(VICONdat{iFile});
        
        % Retrieve analogs
        VICONdat{iFile} = getVICONAnalogIn(VICONdat{iFile}, analogs);

        %Retrieve KUKA position
        %VICONdat{iFile} = getKUKAposition(VICONdat{iFile});

        % Syncronization
        % Number of square peaks in TDT trig in VICON
        %VICONdat{iFile}.TDT_Trigger.times = alignTDT2VICON_robust(TDTdat{iFile}, VICONdat{iFile}, 'Brienne');
         [VICONdat{iFile}] = alignTDT2VICON_allcases(TDTdat{iFile}, VICONdat{iFile}, 'Brienne', TDTdat{iFile}.date);

        [peak, VICONdat{iFile}.BKR_trigger.samples] = (findpeaks(VICONdat{iFile}.analog.data(:,1), 'MinPeakHeight', 1, 'MinPeakProminence', 1));
        VICONdat{iFile}.BKR_trigger.times = VICONdat{iFile}.BKR_trigger.samples(1)/VICONdat{iFile}.analog.framerate;

        % Retrieve events
        VICONdat{iFile} = getVICONEvents(VICONdat{iFile});
        
        
    end
    
    if any(strcmp(recsystems,'BLACKROCK'))
        BKRdat{iFile}.TDT_Trigger.times = VICONdat{iFile}.TDT_Trigger.times - VICONdat{iFile}.BKR_trigger.times;
        BKRdat{iFile}.TDT_Trigger.samples = floor(BKRdat{iFile}.TDT_Trigger.times * BKRdat{iFile}.nsData.MetaTags.SamplingFreq);
        BKRdat{iFile}.FR = preprocessNeuralData(BKRdat{iFile}.nevData, artremovalparams);

    end
end

%% Additional processing (more advanced and optional)

for iFile = 1 : length(TDTdat)
    iFile
     VICONdat{iFile} = getVICONEvents(VICONdat{iFile});
    % Cleaning events
    %[TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile});
    if TDTdat{iFile}.isstim == 1 || TDTdat{iFile}.isstim == 0
        [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile},...
            {'startMov', 'grabMov', 'endMov'});
    else
        [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile},...
            {'startMov', 'grabMov', 'endMov'});
    end
    % Extracting movement attempt indicator
    TDTdat{iFile} = extractMovementAttempt(TDTdat{iFile}, 0.5);
    % Transforming triggers in blackrock in burst triggers
    if TDTdat{iFile}.isstim == 1
        BKRdat{iFile} = extractBursts_BKR(BKRdat{iFile});
        % Clean BKR bursts
        BKRdat{iFile} = cleanBursts_BKR(BKRdat{iFile}, Animal, TDTdat{iFile}.date);
    end
    % Clean force signals 
    VICONdat{iFile} = cleanForceMkBr(VICONdat{iFile});
   
    %Check stim timing
    if TDTdat{iFile}.isstim == 1
        isstimtimed = checkStimTiming(BKRdat{iFile}); 
        TDTdat{iFile}.isstimtimed = isstimtimed;
        stimelecs = checkStimElectrode(BKRdat{iFile}); 
        TDTdat{iFile}.stimelecs  = stimelecs;
    end
    % Assigned success (to do just if you don't filter for success)
    VICONdat{iFile} = assignSuccess(VICONdat{iFile});
end


%% Plot means
nP = 1000;
allFRS1 = zeros(48, nP);
allFRM1 = zeros(48, nP);
allSHO = [];
allELB = [];
allWRI = [];
allWRIx = [];
allELBz = [];
allgrasp = [];
cortexmapping = load('C:\GIT\unifr_cervicalEES\Structs\Ygritte_ArrayMapIndices.mat')
for iFile = 1 : length(TDTdat)
    for iM = 1 : length(VICONdat{iFile}.events_tokeep.startMov)
        eventsIdxsTDT(1) = floor(TDTdat{iFile}.events_tokeep.startMov(iM)*TDTdat{iFile}.fs);
        eventsIdxsTDT(2) = floor(TDTdat{iFile}.events_tokeep.endMov(iM)*TDTdat{iFile}.fs);
        eventsIdxsTDT(3) = floor(TDTdat{iFile}.events_tokeep.grabMov(iM)*TDTdat{iFile}.fs);
        grabTDT = (eventsIdxsTDT(3)  - eventsIdxsTDT(1))/TDTdat{iFile}.fs; 

        % VICON
        eventsIdxsVICON(1) = floor(VICONdat{iFile}.events_tokeep.startMov(iM)*100);
        eventsIdxsVICON(2) = floor(VICONdat{iFile}.events_tokeep.endMov(iM)*100);
        eventsIdxsVICON(3) = floor(VICONdat{iFile}.events_tokeep.grabMov(iM)*100);
        grabVICON = (eventsIdxsVICON(3)  - eventsIdxsVICON(1)); 

        % BKR 
        eventsIdxsBKR(1) = floor(BKRdat{iFile}.events_tokeep.startMov(iM)*100);
        eventsIdxsBKR(2) = floor(BKRdat{iFile}.events_tokeep.endMov(iM)*100);
        eventsIdxsBKR(3) = floor(BKRdat{iFile}.events_tokeep.grabMov(iM)*100);
        grabBKR = eventsIdxsBKR(3)  - eventsIdxsBKR(1); 
        
       
        
        % Firingrates
        FRS1 = BKRdat{iFile}.FR(cortexmapping.s1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 
        
        FRM1 = BKRdat{iFile}.FR(cortexmapping.m1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 
        
        if iFile ==1 && iM == 1
            [sortedFRS1, Is] = sortNeuralChannels(FRS1);
            [sortedFRM1, Im] = sortNeuralChannels(FRM1);
        else
            sortedFRS1 = remapchannels(FRS1,Is); 
            sortedFRM1  = remapchannels(FRM1,Im) ;
        end
        
        
        oldtime = linspace(1, size(sortedFRS1, 2), size(sortedFRS1, 2));
        newtime = linspace(1, size(sortedFRS1, 2), nP);
        for i = 1 : size(sortedFRS1, 1)
            intsortedFRS1(i,:) = interp1(oldtime, sortedFRS1(i,:), newtime);
        end
        for i = 1 : size(sortedFRM1, 1)
            intsortedFRM1(i,:) = interp1(oldtime, sortedFRM1(i,:), newtime);
        end
        allFRS1 = allFRS1+ intsortedFRS1; 
        allFRM1 = allFRM1+ intsortedFRM1; 
        
        % kinematics
        oldtime = linspace(1, length(VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2))),...
            length(VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2))));
        newtime =  linspace(1, length(VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2))),...
            nP);
        
         interpSHO = interp1(oldtime, ...
            VICONdat{iFile}.jointAngles.jointangles.Shoulder(eventsIdxsVICON(1) :eventsIdxsVICON(2)), ...
            newtime);
        interpELB = interp1(oldtime, ...
            VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2)), ...
            newtime);
        interpWRI = interp1(oldtime, ...
            VICONdat{iFile}.jointAngles.jointangles.Wrist(eventsIdxsVICON(1) :eventsIdxsVICON(2)), ...
            newtime);
        
        interpWRIx = interp1(oldtime, ...
            VICONdat{iFile}.jointAngles.wristpos.x(eventsIdxsVICON(1) :eventsIdxsVICON(2)), ...
            newtime);
        interpELBz= interp1(oldtime, ...
            VICONdat{iFile}.kinematic.z(eventsIdxsVICON(1) :eventsIdxsVICON(2), 3)', ...
            newtime);
        allSHO = [allSHO; interpSHO];
        allELB = [allELB; interpELB];
        allWRI = [allWRI; interpWRI];
        allWRIx = [allWRIx; interpWRIx];
        allELBz = [allELBz; interpELBz];
        
        allgrasp = [allgrasp; grabVICON*nP/length(oldtime)];
    end
end
%[allFRS1_s, I] = sortNeuralChannels(allFRS1);
%[allFRM1_s, I] = sortNeuralChannels(allFRM1);



ntrials = size(allELBz, 1);
for i=1:size(allFRM1, 1)
   allFRM1(i,:)=allFRM1(i,:)./ntrials;
end
[max_mat,idx_max]=max(allFRM1,[],2);
for i=1:size(allFRM1, 1)
   allFRM1_norm(i,:)=allFRM1(i,:)./(max_mat(i));
end

figure
subplot(4, 1, 1)
imagesc(allFRS1./ntrials)

subplot(4, 1, 2)
imagesc(allFRM1_norm)
colormap('bone')
subplot(4, 1, 3)
hold on
stdshade(allSHO, 0.2, [1 1 0])
stdshade(allELB, 0.2, [1 0 0])
stdshade(allWRI, 0.2, [0 0 1])
plot([mean(allgrasp), mean(allgrasp)], [50 200], '--k')

subplot(4, 1, 4)
hold on
stdshade(allWRIx, 0.2, [0 0 1])
stdshade(allELBz, 0.2, [1 0 0])

plot([mean(allgrasp), mean(allgrasp)], [50 200], '--k')


%% EMG 
for iFile = 1 : length(TDTdat)
    for iC = 1 : size(TDTdat{iFile}.EMGb, 2)
        %TDTdat{iFile}.envEMGb(:,iC) = envelope(TDTdat{iFile}.EMGb(:,iC), 30,'peak');
        TDTdat{iFile}.envEMGb(:,iC) = highPass(TDTdat{iFile}.EMGb(:,iC), 12207.03, 0.1, 3);
        TDTdat{iFile}.envEMGb(:,iC) = abs(TDTdat{iFile}.envEMGb(:,iC));
        TDTdat{iFile}.envEMGb(:,iC) = lowPass(TDTdat{iFile}.envEMGb(:,iC), 12207.03, 20, 3);
    end
end


nP = 500;
interpEMG= cell(8,1);
for iC = 1 : length(interpEMG)
    interpEMG{iC} = [];
end
for iFile = 1 : length(TDTdat)
    for iM = 1 : length(VICONdat{iFile}.events_tokeep.startMov)
        eventsIdxsTDT(1) = floor(TDTdat{iFile}.events_tokeep.startMov(iM)*TDTdat{iFile}.fs);
        eventsIdxsTDT(2) = floor(TDTdat{iFile}.events_tokeep.endMov(iM)*TDTdat{iFile}.fs);
        eventsIdxsTDT(3) = floor(TDTdat{iFile}.events_tokeep.grabMov(iM)*TDTdat{iFile}.fs);
        grabTDT = (eventsIdxsTDT(3)  - eventsIdxsTDT(1))/TDTdat{iFile}.fs; 

        
        % Envelopes
        oldtime = linspace(1, eventsIdxsTDT(2) - eventsIdxsTDT(1) + 1,...
            eventsIdxsTDT(2) - eventsIdxsTDT(1) + 1 );
        newtime =  linspace(1, eventsIdxsTDT(2) - eventsIdxsTDT(1) + 1,...
            nP);
        for iC = 1 : size(TDTdat{iFile}.EMGb, 2)
            temp = interp1(oldtime, ...
                TDTdat{iFile}.envEMGb(eventsIdxsTDT(1) :eventsIdxsTDT(2),iC), ...
                newtime);
            interpEMG{iC} = [interpEMG{iC}, temp'];
        end
        
    end
end

figure
for iC = 1 : length(interpEMG)
    subplot(8, 1, iC)
    stdshade(interpEMG{iC}', 0.2, [0 0 0])
    title(TDTdat{1}.muscles{iC})
end
%% plot all together for specific trial
cortexmapping = load('C:\GIT\unifr_cervicalEES\Structs\Ygritte_ArrayMapIndices.mat')

iFile = 1;
iM = 3;
% TDT
eventsIdxsTDT(1) = floor(TDTdat{iFile}.events_tokeep.startMov(iM)*TDTdat{iFile}.fs);
eventsIdxsTDT(2) = floor(TDTdat{iFile}.events_tokeep.endMov(iM)*TDTdat{iFile}.fs);
eventsIdxsTDT(3) = floor(TDTdat{iFile}.events_tokeep.grabMov(iM)*TDTdat{iFile}.fs);
grabTDT = (eventsIdxsTDT(3)  - eventsIdxsTDT(1))/TDTdat{iFile}.fs; 

% VICON
eventsIdxsVICON(1) = floor(VICONdat{iFile}.events_tokeep.startMov(iM)*100);
eventsIdxsVICON(2) = floor(VICONdat{iFile}.events_tokeep.endMov(iM)*100);
eventsIdxsVICON(3) = floor(VICONdat{iFile}.events_tokeep.grabMov(iM)*100);
grabVICON = (eventsIdxsVICON(3)  - eventsIdxsVICON(1))/100; 

% BKR 
eventsIdxsBKR(1) = floor(BKRdat{iFile}.events_tokeep.startMov(iM)*100);
eventsIdxsBKR(2) = floor(BKRdat{iFile}.events_tokeep.endMov(iM)*100);
eventsIdxsBKR(3) = floor(BKRdat{iFile}.events_tokeep.grabMov(iM)*100);
grabBKR = eventsIdxsBKR(3)  - eventsIdxsBKR(1); 
% Computing quantities to plot
FRS1 = BKRdat{iFile}.FR(cortexmapping.s1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 
[sortedFRS1, I] = sortNeuralChannels(FRS1);
FRM1 = BKRdat{iFile}.FR(cortexmapping.m1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 
[sortedFRM1, I] = sortNeuralChannels(FRM1);

figure
nplots = 6;
% Plot 1 : cortical activity of S1
subplot(nplots, 2, [1 2])
imagesc(sortedFRS1)
hold on
plot([grabBKR grabBKR], [1 size(sortedFRS1, 1)], '--r')
colorbar

% Plot 2 : cortical activity of M1
subplot(nplots, 2, [3 4])
imagesc(sortedFRM1)
hold on
plot([grabBKR grabBKR], [1 size(sortedFRM1, 1)], '--r')
colorbar
colormap('bone')

% Plot 3 : KINE
% Angles
subplot(nplots, 2, [5 6])
timeKINE = linspace(0,...
    length(VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2)))*1/100, ...
    length(VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2))));
hold on
plot(timeKINE, VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2)), 'r')
plot(timeKINE, VICONdat{iFile}.jointAngles.jointangles.Shoulder(eventsIdxsVICON(1) :eventsIdxsVICON(2)), 'y')
plot(timeKINE, VICONdat{iFile}.jointAngles.jointangles.Wrist(eventsIdxsVICON(1) :eventsIdxsVICON(2)), 'b')
hold on
plot([grabVICON grabVICON], [50 200], '--r')

% Trajectories
subplot(nplots, 2, [7 8])
hold on
plot(timeKINE, ...
    VICONdat{iFile}.jointAngles.wristpos.x(eventsIdxsVICON(1):eventsIdxsVICON(2)) - VICONdat{iFile}.jointAngles.wristpos.x(eventsIdxsVICON(1)), 'b')
plot(timeKINE, VICONdat{iFile}.kinematic.z(eventsIdxsVICON(1) :eventsIdxsVICON(2), 3) - VICONdat{iFile}.kinematic.z(eventsIdxsVICON(1), 3), 'r')
plot([grabVICON grabVICON], [-200 200], '--r')

% Plot 4 : EMG
subplot(nplots, 2, [9 10])
hold on
limsplotEMG = sum(max(TDTdat{iFile}.EMGb))/2;
spacingEMG = limsplotEMG/(size(TDTdat{iFile}.EMGb, 2)) ;
timeEMG = linspace(0,...
    size(TDTdat{iFile}.EMGb(eventsIdxsTDT(1) :eventsIdxsTDT(2) ,1), 1)*1/12207.03, ...
    size(TDTdat{iFile}.EMGb(eventsIdxsTDT(1) :eventsIdxsTDT(2) ,1), 1));
% Plotting EMGs
for iC = 1:size(TDTdat{iFile}.EMGb, 2)-1
    plot(timeEMG, TDTdat{iFile}.EMGb(eventsIdxsTDT(1) :eventsIdxsTDT(2) ,iC) - (iC-1)*spacingEMG, 'k')
end
plot([grabTDT grabTDT], [-limsplotEMG 0], '--r')

%completetrialPlot(TDTdat{iFile},VICONdat{iFile},BKRdat{iFile}, {'EMG'})


clear stickshoulder stickelbow stickwrist
stickshoulder(:,1) = VICONdat{iFile}.kinematic.x(:, 1);
stickshoulder(:,2) = VICONdat{iFile}.kinematic.z(:, 1);
stickelbow(:,1) = VICONdat{iFile}.kinematic.x( :, 3);
stickelbow(:,2) = VICONdat{iFile}.kinematic.z( :,3);
stickwrist(:,1) = VICONdat{iFile}.kinematic.x(:, 5);
stickwrist(:,2) = VICONdat{iFile}.kinematic.z(:, 5);
subplot(nplots, 2, [11])
plotsagittalStickDiagram(stickshoulder, stickelbow, stickwrist, ...
    [eventsIdxsVICON(1) eventsIdxsVICON(3) eventsIdxsVICON(3)], 0, 3, ...
    repmat([1 0 0], 100, 1),  repmat([0 0 1], 100, 1))

subplot(nplots, 2, [12])
plotsagittalStickDiagram(stickshoulder, stickelbow, stickwrist, ...
    [eventsIdxsVICON(3) eventsIdxsVICON(2) eventsIdxsVICON(2)], 0, 3, ...
    repmat([1 0 0], 100, 1),  repmat([0 0 1], 100, 1))


%% just stick
figure

clear stickshoulder stickelbow stickwrist
stickshoulder(:,1) = lowPass(VICONdat{iFile}.kinematic.x(:, 1), 100, 6, 3);
stickshoulder(:,2) = lowPass(VICONdat{iFile}.kinematic.z(:, 1), 100, 6, 3);
stickelbow(:,1) =  lowPass(VICONdat{iFile}.kinematic.x( :, 3), 100, 6, 3);
stickelbow(:,2) =  lowPass(VICONdat{iFile}.kinematic.z( :,3), 100, 6, 3);
stickwrist(:,1) =  lowPass(VICONdat{iFile}.kinematic.x(:, 5), 100, 6, 3);
stickwrist(:,2) =  lowPass(VICONdat{iFile}.kinematic.z(:, 5), 100, 6, 3);
subplot(1, 2, [1])
plotsagittalStickDiagram(stickshoulder, stickelbow, stickwrist, ...
    [eventsIdxsVICON(1) eventsIdxsVICON(3) eventsIdxsVICON(3)], 0, 2, ...
    repmat([1 0 0], 100, 1),  repmat([0 0 1], 100, 1))
axis equal
subplot(1, 2, [2])
plotsagittalStickDiagram(stickshoulder, stickelbow, stickwrist, ...
    [eventsIdxsVICON(3) eventsIdxsVICON(2) eventsIdxsVICON(2)], 0, 2, ...
    repmat([1 0 0], 100, 1),  repmat([0 0 1], 100, 1))
axis equal


%%
iF = 1; 
iM = 1; 
eventsIdxsTDT(1) = floor(TDTdat{iFile}.events_tokeep.startMov(iM)*TDTdat{iFile}.fs);
eventsIdxsTDT(2) = floor(TDTdat{iFile}.events_tokeep.endMov(iM)*TDTdat{iFile}.fs);
eventsIdxsTDT(3) = floor(TDTdat{iFile}.events_tokeep.grabMov(iM)*TDTdat{iFile}.fs);
grabTDT = (eventsIdxsTDT(3)  - eventsIdxsTDT(1))/TDTdat{iFile}.fs; 

% VICON
eventsIdxsVICON(1) = floor(VICONdat{iFile}.events_tokeep.startMov(iM)*100);
eventsIdxsVICON(2) = floor(VICONdat{iFile}.events_tokeep.endMov(iM)*100);
eventsIdxsVICON(3) = floor(VICONdat{iFile}.events_tokeep.grabMov(iM)*100);
grabVICON = (eventsIdxsVICON(3)  - eventsIdxsVICON(1)); 

% BKR 
eventsIdxsBKR(1) = floor(BKRdat{iFile}.events_tokeep.startMov(iM)*100);
eventsIdxsBKR(2) = floor(BKRdat{iFile}.events_tokeep.endMov(iM)*100);
eventsIdxsBKR(3) = floor(BKRdat{iFile}.events_tokeep.grabMov(iM)*100);
grabBKR = eventsIdxsBKR(3)  - eventsIdxsBKR(1); 



% Firingrates
FRS1 = BKRdat{iFile}.FR(cortexmapping.s1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 

FRM1 = BKRdat{iFile}.FR(cortexmapping.m1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 
[FSM1_s, I] = sortNeuralChannels(FRM1)
figure
subplot(1, 3, 1)
imagesc(FRM1)
subplot(1, 3, 2)
imagesc(FSM1_s)
subplot(1, 3, 3)
sortedmat = remapchannels(FRM1,I) 
imagesc(sortedmat)

%%
% Computing quantities to plot
FRS1 = BKRdat{iFile}.FR(cortexmapping.s1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 
[sortedFRS1, I] = sortNeuralChannels(FRS1);
FRM1 = BKRdat{iFile}.FR(cortexmapping.m1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 
[sortedFRM1, I] = sortNeuralChannels(FRM1);

figure
nplots = 5;
% Plot 1 : cortical activity of S1
subplot(nplots, 1, 1)
imagesc(sortedFRS1)
hold on
plot([grabBKR grabBKR], [1 size(sortedFRS1, 1)], '--r')
colorbar

% Plot 2 : cortical activity of M1
subplot(nplots, 1, 2)
imagesc(sortedFRM1)
hold on
plot([grabBKR grabBKR], [1 size(sortedFRM1, 1)], '--r')
colorbar
colormap('bone')

% Plot 3 : KINE
% Angles
subplot(nplots, 1, 3)
timeKINE = linspace(0,...
    length(VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2)))*1/100, ...
    length(VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2))));
hold on
plot(timeKINE, VICONdat{iFile}.jointAngles.jointangles.Elbow(eventsIdxsVICON(1) :eventsIdxsVICON(2)), 'r')
plot(timeKINE, VICONdat{iFile}.jointAngles.jointangles.Shoulder(eventsIdxsVICON(1) :eventsIdxsVICON(2)), 'y')
plot(timeKINE, VICONdat{iFile}.jointAngles.jointangles.Wrist(eventsIdxsVICON(1) :eventsIdxsVICON(2)), 'b')
hold on
plot([grabVICON grabVICON], [50 200], '--r')

% Trajectories
subplot(nplots, 1, 4)
hold on
plot(timeKINE, ...
    VICONdat{iFile}.jointAngles.wristpos.x(eventsIdxsVICON(1):eventsIdxsVICON(2)) - VICONdat{iFile}.jointAngles.wristpos.x(eventsIdxsVICON(1)), 'b')
plot(timeKINE, VICONdat{iFile}.kinematic.z(eventsIdxsVICON(1) :eventsIdxsVICON(2), 3) - VICONdat{iFile}.kinematic.z(eventsIdxsVICON(1), 3), 'r')
plot([grabVICON grabVICON], [-200 200], '--r')

% Plot 4 : EMG
subplot(nplots, 1, 5)
hold on
limsplotEMG = sum(max(TDTdat{iFile}.EMGb))/2;
spacingEMG = limsplotEMG/(size(TDTdat{iFile}.EMGb, 2)) ;
timeEMG = linspace(0,...
    size(TDTdat{iFile}.EMGb(eventsIdxsTDT(1) :eventsIdxsTDT(2) ,1), 1)*1/12207.03, ...
    size(TDTdat{iFile}.EMGb(eventsIdxsTDT(1) :eventsIdxsTDT(2) ,1), 1));
% Plotting EMGs
for iC = 1:size(TDTdat{iFile}.EMGb, 2)-1
    plot(timeEMG, TDTdat{iFile}.EMGb(eventsIdxsTDT(1) :eventsIdxsTDT(2) ,iC) - (iC-1)*spacingEMG, 'k')
end
plot([grabTDT grabTDT], [-limsplotEMG 0], '--r')

%%
shoulder.x = VICONdat{iFile}.kinematic.x(:,1);
shoulder.y = VICONdat{iFile}.kinematic.z(:,1);

elbow.x = VICONdat{iFile}.kinematic.x(:,2);
elbow.y = VICONdat{iFile}.kinematic.z(:,2);

wrist.x = VICONdat{iFile}.kinematic.x(:,3);
wrist.y = VICONdat{iFile}.kinematic.z(:,3);

shift = 0.02;
Dt = 2;
iFile = 1;
iM = 2;
figure
hold on
eventsIdxs(1) = floor(VICONdat{1}.events_tokeep.startMov(iM)*100);
eventsIdxs(2) = floor(VICONdat{1}.events_tokeep.endMov(iM)*100);
stickdiagram_ygritte_behav(shoulder, elbow, wrist, eventsIdxs, shift, Dt)
%%
cortexmapping = load('C:\GIT\unifr_cervicalEES\Structs\Ygritte_ArrayMapIndices.mat')
iFile = 4;
iM = 1;
eventsIdxsBKR(1) = floor(BKRdat{iFile}.events_tokeep.startMov(iM)*100);
eventsIdxsBKR(2) = floor(BKRdat{iFile}.events_tokeep.endMov(iM)*100);
grab = floor(BKRdat{iFile}.events_tokeep.grabMov(iM)*100) - eventsIdxsBKR(1);
figure
FRS1 = BKRdat{iFile}.FR(cortexmapping.s1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 
[sortedFRS1, I] = sortNeuralChannels(FRS1);
FRM1 = BKRdat{iFile}.FR(cortexmapping.m1_idx, eventsIdxsBKR(1): eventsIdxsBKR(2)); 
[sortedFRM1, I] = sortNeuralChannels(FRM1);
figure
subplot(2, 1, 1)
imagesc(sortedFRS1)
hold on
plot([grab grab], [1 size(sortedFRS1, 1)], '--r')
colorbar
subplot(2, 1, 2)
imagesc(sortedFRM1)
hold on
plot([grab grab], [1 size(sortedFRS1, 1)], '--r')
colormap('bone')
colorbar


%% Performance in baseline
nB = 0;
nS = 0;
for iFile = 1 : length(VICONdat)
    temp = length(find(cellfun(@(x) strcmp(x, 'BadTrial'), VICONdat{iFile}.event.labels)));
    nB = nB + temp;
    
    temp = length(VICONdat{iFile}.events_tokeep.startMov);
    nS = nS + temp;
end
tot = nS + nB;
perS = nS/tot*100 
perB = nB/tot*100

figure
pie([perS, perB])