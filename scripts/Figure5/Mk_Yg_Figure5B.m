warning('off','all')
clear all
clc
close all
Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES';
load('/Users/barrab01/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')

Animal = 'Ygritte';
%expDates = {'20190608', '20190611',  '20190614', '20190624'};
%expDates = { '20190802', '20190809', '20190812', '20190813'}; % 
% selectedFiles = {
%     [1 2]; ...
%     [ 4 5 6 7 10 14 15 16 ];...
%     [ 2 8 10  12 13  17 18]; ...% add the 14
%     [2 3 4 10 11 12 13 14 15 16]; ...%2    
% }% each cell is a date
% stimMask = {
%     [2 2]; ...
%     [0 0 1 1 1 1 1 1]; ...
%     [1 1 1 1  1  1 1]; ...
%     [0 0 1 1 1 1 1 1 1 1]; ...%0  
% }% each cell is a date, 1 means stim file, 0 is baseline

% ---------------------------------------
% expDates = {'20190802', '20190813'};
% selectedFiles = {
%     [1 2]; ...
%     [2 3 4  10 11 12 13 14 15 16  ]; ...%2    
% }% each cell is a date
% stimMask = {
%      [2 2]; ...
%     [ 0 0 0   1 1 1 1 1 1 1  ]; ...%0  
% }% each cell is a date, 1 means 
% ---------------------------------------
% expDates = {'20190802', '20190809'};
% selectedFiles = {
%     [1 2]; ...
%     [4  6 7 10  12 14 15 16 17]; ...%2    
% }% each cell is a date
% stimMask = {
%      [2 2]; ...
%     [ 0 1 1  1 1 1 1 1 1 ]; ...%0  
% }% each cell is a date, 1 means 
% ---------------------------------------


% expDates = {'20190802', '20190812'};
% selectedFiles = {
%     [1 2]; ...
%     [2  5 12 13   17 18 ]; ...%2    
% }% each cell is a date
% stimMask = {
%      [2 2]; ...
%     [ 0  1 1 1  1 1]; ...%0  
% }% each cell is a date, 1 means 
% ---------------------------------------

expDates = {'20190802', '20190820'};
selectedFiles = {
    [1 2]; ...
    [1 3 4 6 8 9 10 11 12 14 17]; ...%2    
}% each cell is a date
stimMask = {
     [2 2]; ...
    [0 0 1 1 1 1 1 1 1 1 1]; ...%0  
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
        if strcmp(expDate, '20190529')
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
        %VICONdat{iFile}.TDT_Trigger.times = alignTDT2VICON_robust(TDTdat{iFile}, VICONdat{iFile}, 'Brienne');
        [VICONdat{iFile}] = alignTDT2VICON_allcases(TDTdat{iFile}, VICONdat{iFile}, 'Ygritte', TDTdat{iFile}.date);

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
     VICONdat{iFile} = getVICONEvents(VICONdat{iFile});
    % Cleaning events
    %[TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile});
    if TDTdat{iFile}.isstim == 1 || TDTdat{iFile}.isstim == 0
        [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile},...
            {'startMov', 'grabMov', 'endMov' , 'reach'});
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

%% Plot of full trial
% effstim = 0;
% for iFile = 5 : length(TDTdat)
% figure
% 
% completetrialPlot(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile})
% TDTdat{iFile}.isstimtimed.reach
% effstim = effstim + length(find(TDTdat{iFile}.isstimtimed.reach));
% pause
% end

%% counting reach
countbasreach= 0;
countstimreach = 0;
countlesreach = 0;
for iFile = 1 : length(TDTdat)

% Check stim timing
    if TDTdat{iFile}.isstim == 0
        countlesreach = countlesreach + length(VICONdat{iFile}.events_tokeep.startMov);
        
    end
    if TDTdat{iFile}.isstim == 1
        countstimreach = countstimreach + length(VICONdat{iFile}.reach);
        
    end
    if TDTdat{iFile}.isstim == 2
        countbasreach = countbasreach + length(VICONdat{iFile}.events_tokeep.grabMov);
        
    end
end
countbasreach
countlesreach
countstimreach
%% Kinematic features computation
        featureslabels = {...
                 'maxSHO', ...% max shoulder angle
                 'maxELB', ...% max elbow angles
                 'minSHO', ...
                 'minELB', ...
                 'excSHO', ...
                 'excELB', ...
                 'maxvelSHO', ...
                 'maxvelELB', ...
                 'avgvelSHO', ...
                 'avgvelELB', ...
                 'maxWRIPx', ...
                 'maxWRIPz', ...
                 'minWRIPx', ...
                 'minWRIPz', ...
                 'maxELBPx', ...
                 'maxELBPz', ...
                 'minELBPx', ...
                 'minELBPz', ...
                 'avgELBPz', ...
                 'maxSHOPx', ...
                 'maxSHOPz', ...
                 'minSHOPx', ...
                 'minSHOPz', ...
                 'maxWRIVx', ...
                 'maxWRIVz', ...
                 'avgWRIVx', ...
                 'avgWRIVz', ...
                 'smoothness', ...
                 'timeREACH', ...%
                 'lengthREACH', ...%             
                 };        
fsa = 1000;
fsk = 100;
features = [];
isstim = [];
isstimtimedreach = [];
isstimtimedpull = [];
issuccess = [];
stimelecs = [];
snips = {};
stimtimes = {};
pairedEMGvalues = [];
elbowheight = [];
for iFile = 1 : length(VICONdat)
    filteredShoulder = lowPass(VICONdat{iFile}.jointAngles.jointangles.Shoulder', ...
        100, 6, 3);
    filteredElbow = lowPass(VICONdat{iFile}.jointAngles.jointangles.Elbow', ...
        100, 6, 3);
    
    
    filteredElbowx = lowPass(VICONdat{iFile}.kinematic.x(:,2), ...
        100, 6, 3);
    filteredElbowy = lowPass(VICONdat{iFile}.kinematic.z(:,2), ...
        100, 6, 3);
    filteredShoulderx = lowPass(-VICONdat{iFile}.kinematic.x(:,1), ...
        100, 6, 3);
    filteredShouldery = lowPass(-VICONdat{iFile}.kinematic.z(:,1), ...
        100, 6, 3);
    filteredwristx = lowPass(VICONdat{iFile}.jointAngles.wristpos.x, ...
        100, 6, 3);
    filteredwristz = lowPass(VICONdat{iFile}.jointAngles.wristpos.z, ...
        100, 6, 3);
    
    filteredElbowHeight = lowPass(VICONdat{iFile}.kinematic.z(:,2), ...
        100, 6, 3);
    
    for iM =1 : length(VICONdat{iFile}.events_tokeep.startMov)
        
        % Start : present everywhere
        startM = VICONdat{iFile}.events_tokeep.startMov(iM);
        if startM == 0
            startM =0.01;
        end
        startM_ia = floor(startM*fsa);
        startM_ik = floor(startM*fsk);
        % End : present everywhere 
        endM = VICONdat{iFile}.events_tokeep.endMov(iM);
        endM_ia = floor(endM*fsa);
        endM_ik = floor(endM*fsk);
        % Grab
        grabM = VICONdat{iFile}.events_tokeep.grabMov(iM);
        grabM_ia = floor(grabM*fsa);
        grabM_ik = floor(grabM*fsk);
        mystart_k =startM_ik;  myend_k =grabM_ik; mygrab_k = grabM_ik;  myend_ksnip =endM_ik;
        mystart_a =startM_ia;  myend_a =grabM_ia; mygrab_a = grabM_ia;
        timeREACH = grabM - startM;
        
        
        % hacking the starts
        
        xwrist = VICONdat{iFile}.jointAngles.wristpos.x(mystart_k:myend_k) - VICONdat{iFile}.jointAngles.wristpos.x(mystart_k);
        vwrist = diff(xwrist);
        [ startidx]= min(find(xwrist<-50));
%         plot(xwrist, 'k')
%         hold on
%         plot(startidx, xwrist(startidx), 'or')
%         pause
%         clf
        mystart_k = mystart_k + startidx;
        startM = mystart_k/VICONdat{iFile}.kinematic.framerate; 
        % Angles
        maxSHO = max(filteredShoulder(mystart_k:myend_k));
        maxELB = max(filteredElbow(mystart_k:myend_k));
        minSHO = min(filteredShoulder(mystart_k:myend_k));
        minELB = min(filteredElbow(mystart_k:myend_k));
        excSHO = maxSHO - minSHO;
        excELB = maxELB - minELB;
        velSHO = diff(filteredShoulder(mystart_k:myend_k));
        velELB = diff(filteredElbow(mystart_k:myend_k));
        maxvelSHO = max(velSHO);
        maxvelELB = max(velELB);
        avgvelSHO = mean(velSHO);
        avgvelELB = mean(velELB);
        %
        maxELBPx = max(filteredElbowx(mystart_k:myend_k) - filteredElbowx(mystart_k));
        maxELBPz = max(filteredElbowy(mystart_k:myend_k) - filteredElbowy(mystart_k));
        minELBPx = min(filteredElbowx(mystart_k:myend_k) - filteredElbowx(mystart_k));
        minELBPz = min(filteredElbowy(mystart_k:myend_k) - filteredElbowy(mystart_k));
        avgELBPz = mean(filteredElbowy(myend_k-10 :myend_k))-filteredElbowy(mystart_k);
        maxSHOPx = max(filteredShoulderx(mystart_k:myend_k) - filteredShoulderx(mystart_k));
        maxSHOPz = max(filteredShouldery(mystart_k:myend_k) - filteredShouldery(mystart_k));
        minSHOPx = min(filteredShoulderx(mystart_k:myend_k) - filteredShoulderx(mystart_k));
        minSHOPz = min(filteredShouldery(mystart_k:myend_k) - filteredShouldery(mystart_k));
        
        
        % Wrist trajectory
        maxWRIPx = max(VICONdat{iFile}.jointAngles.wristpos.x(mystart_k:myend_k)...
            - VICONdat{iFile}.jointAngles.wristpos.x(mystart_k));
        maxWRIPz = max(VICONdat{iFile}.jointAngles.wristpos.z(mystart_k:myend_k)-...
            VICONdat{iFile}.jointAngles.wristpos.z(mystart_k));
        minWRIPx = min(VICONdat{iFile}.jointAngles.wristpos.x(mystart_k:myend_k)- ...
            VICONdat{iFile}.jointAngles.wristpos.x(mystart_k));
        minWRIPz = min(VICONdat{iFile}.jointAngles.wristpos.z(mystart_k:myend_k)- ...
            VICONdat{iFile}.jointAngles.wristpos.z(mystart_k));
        
        WRIVx = diff(VICONdat{iFile}.jointAngles.wristpos.x(mystart_k:myend_k));
        WRIVz = diff(VICONdat{iFile}.jointAngles.wristpos.z(mystart_k:myend_k));
        
        maxWRIVx = max(WRIVx);
        maxWRIVz = max(WRIVz);
        avgWRIVx = mean(WRIVx);
        avgWRIVz = mean(WRIVz);
        
        % Smoothness
        xwrist = VICONdat{iFile}.jointAngles.wristpos.x(mystart_k:myend_k);
        zwrist = VICONdat{iFile}.jointAngles.wristpos.z(mystart_k:myend_k);
        lengthMOV = sum(sqrt(diff(xwrist).^2 + diff(zwrist).^2 ));
        jerkx = diff(diff(diff(xwrist))); jerkz = diff(diff(diff(zwrist))); 
        jerk = sqrt(jerkx.^2+ jerkz.^2);
        durMOV = grabM - startM;
        smoothness = sqrt(1/3*sum(abs(jerk))*(durMOV.^5/lengthMOV.^2));
       
        
        % Movement trajectories length
        xwrist = VICONdat{iFile}.jointAngles.wristpos.x(mystart_k:mygrab_k);
        zwrist = VICONdat{iFile}.jointAngles.wristpos.z(mystart_k:mygrab_k);
        lengthREACH = sum(sqrt(diff(xwrist).^2 + diff(zwrist).^2 ));
       
       featvect = [maxSHO, ...% max shoulder angle
         maxELB, ...% max elbow angles
         minSHO, ...
         minELB, ...
         excSHO, ...
         excELB, ...
         maxvelSHO, ...
         maxvelELB, ...
         avgvelSHO, ...
         avgvelELB, ...
         maxWRIPx, ...
         maxWRIPz, ...
         minWRIPx, ...
         minWRIPz, ...
         maxELBPx, ...
         maxELBPz, ...
         minELBPx, ...
         minELBPz, ...
         avgELBPz, ...
         maxSHOPx, ...
         maxSHOPz, ...
         minSHOPx, ...
         minSHOPz, ...
         maxWRIVx, ...
         maxWRIVz, ...
         avgWRIVx, ...
         avgWRIVz, ...
         smoothness, ...
         timeREACH, ...
         lengthREACH, ...%
         ];
       timestartstim = [];
        % Finding stim time
        if TDTdat{iFile}.isstim ==1
            for iE = 3% 1 : length(BKRdat{iFile}.burstsStart)
                allstarts = BKRdat{iFile}.burstsStart{iE} + VICONdat{iFile}.BKR_trigger.times;
                
                idx_burst = intersect(find(allstarts>startM), find(allstarts<grabM));
                timestartstim = [timestartstim; allstarts(idx_burst) - startM];
                
            end
        end
       
        %Possibilities for the snips
        % WRIST x
        myend_ksnip = myend_k;
        snip = filteredwristx(mystart_k:myend_ksnip) - filteredwristx(mystart_k);
        % WRIST z
        snip = filteredwristz(mystart_k:myend_ksnip) - filteredwristz(mystart_k);
        % ELBOW angle  
        snip =filteredElbow(mystart_k:myend_ksnip) - filteredElbow(mystart_k);
        % Shoulder x
        %snip =filteredShoulder(mystart_k:myend_ksnip) - filteredShoulder(mystart_k);
        % Elbow height
         snip = filteredElbowHeight(mystart_k:myend_ksnip) - filteredElbowHeight(mystart_k);
        elbowheight = [elbowheight; max(snip)]; 
        %snip = diff(snip); 
          
        sniptimelength = length(snip)*(1/VICONdat{1}.kinematic.framerate);
        
        t1 = linspace(1, length(snip), length(snip));
        t2 =  linspace(1, length(snip),1000);
        snip = interp1(t1, snip, t2);
        if  timestartstim/sniptimelength >1
            disp('problema')
            pause
        end
        
        samplestartstim = floor(1000*timestartstim/sniptimelength);
        %snip = diff(snip);
        if  TDTdat{iFile}.isstim == 1
            if  TDTdat{iFile}.isstimtimed.reach(iM) == 1
                % find the moment the stim intervenes
                startMTDT = TDTdat{iFile}.events_tokeep.startMov(iM); 
                samplestartMTDT =floor( startMTDT*TDTdat{iFile}.fs);
                grabMTDT = TDTdat{iFile}.events_tokeep.grabMov(iM); 
                samplegrabMTDT =floor( grabMTDT*TDTdat{iFile}.fs);
                timestartstimTDT = (startM +timestartstim) - VICONdat{iFile}.TDT_Trigger.times;
                samplestartstimTDT = floor( timestartstimTDT*TDTdat{iFile}.fs);
                idxmuscle = 3; 
                TRInostim = (TDTdat{iFile}.EMGb(samplestartMTDT:samplestartstimTDT,idxmuscle)) 
                energyTRInostim = sum(TRInostim.^2)/length(TRInostim);
                TRIstim = TDTdat{iFile}.EMGb(samplestartstimTDT:samplegrabMTDT,idxmuscle);
                energyTRIstim = sum(TRIstim.^2)/length(TRIstim);
                pairedEMGvalues = [pairedEMGvalues;energyTRInostim, energyTRIstim ]; 
            end
        end
         
         if (isreal(featvect))
             disp([num2str(iFile),', ', num2str(iM), ' is real']);
             
             if TDTdat{iFile}.isstim == 1
                 if TDTdat{iFile}.isstimtimed.pull(iM) == 0 &&  TDTdat{iFile}.isstimtimed.reach(iM) == 0
                     continue;
                 end
             end
             features = [features; featvect];
             snips{end + 1} = snip;
             stimtimes{end+1} = samplestartstim;
            if samplestartstim > 1000
                disp('problema anche qui')
                pause
            end
             isstim  = [isstim; TDTdat{iFile}.isstim];
             if (VICONdat{iFile}.assignedSuccess(iM) == 0)
                 issuccess = [issuccess; 0];
             else
                 issuccess = [issuccess; 1];
             end
                 
             try
                isstimtimedpull  = [isstimtimedpull; TDTdat{iFile}.isstimtimed.pull(iM)];
             catch
                 isstimtimedpull  = [isstimtimedpull; 0];
             end
             try
                isstimtimedreach  = [isstimtimedreach; TDTdat{iFile}.isstimtimed.reach(iM)];
             catch
                 isstimtimedreach  = [isstimtimedreach; 0];
             end
             try 
                 stimelecs = [stimelecs;  TDTdat{iFile}.stimelecs]; 
             catch
                 stimelecs = [stimelecs;  zeros(1, 6)]; 
             end
         else
             disp([num2str(iFile),', ', num2str(iM), ' is complex']);
             features = [features; zeros(size(featvect))];
             isstim  = [isstim; NaN];
             isstimtimed  = [isstimtimed; NaN];
         end
    end  
end

i_eliminate = find(isnan(isstim));
features(i_eliminate, :) = [];
isstim(i_eliminate, :) = [];
isstimtimed(i_eliminate, :) = [];


% Snips shaded area
% FORCE x and 3d force no good
% WRIST x is good
BASidx = find(isstim == 2);
LESIONidx = find(isstim == 0);
STIMidx = intersect(find(isstimtimedreach==1), find(stimelecs(:,3) == 1));



%%
nFeat = size(features, 2);
dims = numSubplots(nFeat);

figure
% testing baseline vas lesion vs stim
for iF = 1 : size(features, 2)
    baseline = features(BASidx, iF);
    lesion = features(LESIONidx, iF);
    stim = features(STIMidx, iF);
    
    data = [baseline; lesion; stim];
    datalabels = [ones(size(baseline)); ...
        2*ones(size(lesion)); ...
        3*ones(size(stim))];
%     [P,ANOVATAB,STATS] = kruskalwallis(data, datalabels);
%     multcompare(STATS)
    p(1) = ranksum(baseline, lesion);
    p(2) = ranksum (lesion, stim);
    p(3) = ranksum (baseline, stim);
    
    subplot(dims(1), dims(2), iF)
    hold on
    b = bar([1 2 3], [mean(baseline), mean(lesion), mean(stim)]);
    er = errorbar([1], mean(baseline), std(baseline), std(baseline));
    er = errorbar([2], mean(lesion), std(lesion), std(lesion));
    er = errorbar([3], mean(stim), std(stim), std(stim));
%     text(1.5, mean([mean(baseline), mean(lesion)]), num2str(p(1)))
%     text(2.5, mean([mean(stim), mean(lesion)]), num2str(p(2)))
    if p(1)<0.05
        text(1.5, mean([mean(baseline), mean(lesion)]), '*')
    end
    if p(2)<0.05
        text(2.5, mean([mean(stim), mean(lesion)]), '*')
    end
    title(featureslabels{iF})
end
%% Selecting a special feature

% Elbow height 
iF = 19; 

%iF = nFeat % ELBpy
%iF = 29 % Fy
figure
baseline = features(BASidx, iF);
    lesion = features(LESIONidx, iF);
    stim = features(STIMidx, iF);
    
    data = [baseline; lesion; stim];
    datalabels = [ones(size(baseline)); ...
        2*ones(size(lesion)); ...
        3*ones(size(stim))];
    figure
    [P,ANOVATAB,STATS] = kruskalwallis(data, datalabels);
    % multcompare(STATS)
    p(1) = ranksum(baseline, lesion)
    p(2) = ranksum (lesion, stim);
    p(3) = ranksum (baseline, stim);
    figure
    hold on
    baseline = baseline(not(isoutlier(baseline)));
    lesion = lesion(not(isoutlier(lesion)));
    stim = stim(not(isoutlier(stim)));
    b = bar([1 2 3], [mean(baseline), mean(lesion), mean(stim)]);
    er = errorbar([1], mean(baseline), 0, std(baseline));%/sqrt(length(baseline)));
    er = errorbar([2], mean(lesion), 0, std(lesion));%/sqrt(length(lesion)));
    er = errorbar([3], mean(stim), 0, std(stim));%/sqrt(length(stim)));
    widthdots =0.5;
    randvector = (rand(size(baseline))-0.5)*widthdots;
    plot(ones(size(randvector)) + randvector, baseline, 'ok')
    randvector = (rand(size(lesion))-0.5)*widthdots;
    plot(2*ones(size(randvector)) + randvector, lesion , 'ok')
    randvector = (rand(size(stim))-0.5)*widthdots;
    plot(3*ones(size(randvector)) + randvector, stim , 'ok')
%     text(1.5, mean([mean(baseline), mean(lesion)]), num2str(p(1)))
%     text(2.5, mean([mean(stim), mean(lesion)]), num2str(p(2)))
    if p(1)<0.05
        text(1.5, mean([mean(baseline), mean(lesion)]), '*')
    end
    if p(2)<0.05
        text(2.5, mean([mean(stim), mean(lesion)]), '*')
    end
    if p(3)<0.05
        text(2, mean([mean(stim), mean(baseline)]), '*')
    end
    title(featureslabels{iF})

    
%% Violin plot

iF = 19;  % Elbow height
%iF = nFeat % Movement length


baseline = features(BASidx, iF);
lesion = features(LESIONidx, iF);
stim = features(STIMidx, iF);

data = [baseline; lesion; stim];
datalabels = [ones(size(baseline)); ...
    2*ones(size(lesion)); ...
    3*ones(size(stim))];
figure
[P,ANOVATAB,STATS] = kruskalwallis(data, datalabels);
% multcompare(STATS)
p(1) = ranksum(baseline, lesion)
p(2) = ranksum (lesion, stim);
p(3) = ranksum (baseline, stim);
figure
hold on
baseline = baseline(not(isoutlier(baseline)));
lesion = lesion(not(isoutlier(lesion)));
stim = stim(not(isoutlier(stim)));

data ={baseline, lesion, stim};
violin(data)

widthdots =0.3;
randvector = (rand(size(baseline))-0.5)*widthdots;
plot(ones(size(randvector)) + randvector, baseline, 'ok')
randvector = (rand(size(lesion))-0.5)*widthdots;
plot(2*ones(size(randvector)) + randvector, lesion , 'ok')
randvector = (rand(size(stim))-0.5)*widthdots;
plot(3*ones(size(randvector)) + randvector, stim , 'ok')
%     text(1.5, mean([mean(baseline), mean(lesion)]), num2str(p(1)))
%     text(2.5, mean([mean(stim), mean(lesion)]), num2str(p(2)))
if p(1)<0.05
    text(1.5, mean([mean(baseline), mean(lesion)]), '*')
end
if p(2)<0.05
    text(2.5, mean([mean(stim), mean(lesion)]), '*')
end
if p(3)<0.05
    text(2, mean([mean(stim), mean(baseline)]), '*')
end
title(featureslabels{iF})
%% PCA

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(zscore(features));


figure('units','normalized','outerposition',[0.3 0.2 0.6 0.8])
subplot(2, 2, 1)
imagesc(COEFF)
set(gca, 'ytick', [1: size(features, 2)], 'yticklabel', featureslabels)

weightsCM = [linspace(13 , 125, 100)', ...
    linspace(46, 125, 100)', ...
    linspace(151, 125, 100)'];

weightsCM = [weightsCM ; [linspace(125, 255, 100)', ...
    linspace(125, 183, 100)', ...
    linspace(125,27, 100)']];
hold on
colormap(weightsCM./255)
%plot(EXPLAINED)
[per, ipc] = (max(-diff(EXPLAINED)))
ipc = 1;
subplot(2 ,2, 2)
hold on
scspapercolors = [230 230 230; 125 125 125; 255 183 27];
countbastrials = 0;
intactcolor = [125 125 125]./255;
lesioncolor= [0 0 0];
stimcolor=[255 183 27]./255;
for ip = 1 : size(SCORE, 1)
    if issuccess(ip) == 1
        marker = 'o';
    else
        marker = 'o';
    end
    if isstim(ip) == 1
            color = stimcolor % YELLOW IS STIM
    elseif isstim(ip) == 2 % prelesion
        
         color = intactcolor; % GREY IS BASELINE PRE LESION
    else 
        color = lesioncolor;% BLACK IS LESION NO STIM
    end
    plot3(SCORE(ip,ipc), SCORE(ip,ipc+1) ,SCORE(ip,ipc+2), marker,...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color, ...
        'Markersize',6)

end

% Compute centroids 


C_BAS = [mean(SCORE(BASidx, ipc)), mean(SCORE(BASidx,ipc+1))] ;
C_STIM = [mean(SCORE(STIMidx, ipc)), mean(SCORE(STIMidx, ipc+1))] ;
C_LESION = [mean(SCORE(LESIONidx, ipc)), mean(SCORE(LESIONidx, ipc+1))] ;
hold on
plot(C_BAS(1), C_BAS(2), 'o', 'MarkerFaceColor', intactcolor,'MarkerEdgeColor', [1 1 1], 'Markersize', 12)
plot(C_STIM(1), C_STIM(2),'o', 'MarkerFaceColor', stimcolor,'MarkerEdgeColor', [255 183 27]./255, 'Markersize', 12)
plot(C_LESION(1), C_LESION(2),'o', 'MarkerFaceColor', lesioncolor,'MarkerEdgeColor',[125 125 125]./255, 'Markersize', 10)
xlabel(['PC1'])
ylabel(['PC2'])
zlabel(['PC3'])


% Euclidean distances
% Compute idx 
clear featdistance_STIMLESION featdistance_BASSTIM featdistance_BASLESION
subplot(2, 2, 4)
hold on

featuresfordistance = zscore(features);
for ibas = 1 : length(BASidx)
    
v1 = featuresfordistance(BASidx(ibas), :);
v2 = mean(featuresfordistance(STIMidx, :));
featdistance_BASSTIM(ibas) = sqrt(sum((v1 - v2).^2)); 
    
end

for ibas = 1 : length(BASidx)
    
    v1 = featuresfordistance(BASidx(ibas), :);
    v2 = mean(featuresfordistance(LESIONidx, :));
    featdistance_BASLESION(ibas) = sqrt(sum((v1 - v2).^2)); 
    
end




randdistr = 0.3*rand(length(featdistance_BASSTIM(:)), 1) - 0.3/2;
plot(ones(size(featdistance_BASSTIM(:)))+randdistr , featdistance_BASSTIM(:), 'o', ...
    'MarkerFaceColor', stimcolor, ...
    'MarkerEdgeColor', stimcolor,... 
    'Markersize', 6)
plot(1 , mean(featdistance_BASSTIM(:)), 'o', ...
    'MarkerFaceColor', [0 0 0], ...
    'MarkerEdgeColor', [0 0 0],... 
    'Markersize', 10)

randdistr = 0.3*rand(length(featdistance_BASLESION(:)), 1) - 0.3/2;
plot(2*ones(size(featdistance_BASLESION(:))) + randdistr, featdistance_BASLESION(:), 'o', ...
    'MarkerFaceColor', lesioncolor, ...
    'MarkerEdgeColor', lesioncolor,... 
    'Markersize', 6)
plot(2 , mean(featdistance_BASLESION(:)), 'o', ...
    'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [1 1 1],... 
    'Markersize', 10)

boxplot([ featdistance_BASSTIM(:), featdistance_BASLESION(:),],  'notch'    ,     'on' )


p = ranksum(featdistance_BASSTIM(:), featdistance_BASLESION(:)); 
title({'Dist of feats from centroids (features zscored)',  ['G1 vs G2 p = ',num2str(p)]} )
ylabel('Euclidean distance (feature space)')
set(gca, 'xtick', [1 2 3], 'xticklabel', {'PRELESION vs STIM', 'PRELESION vs NOSTIM', 'STIM vs NOSTIM'})
xlim([0.6 3.2])




