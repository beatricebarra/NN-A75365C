warning('off','all')
clear all
clc
close all
Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES';
load('/Users/barrab01/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')


Animal = 'Brienne';
%expDates = {'20190608', '20190611',  '20190614', '20190624'};
expDates = { '20190529', '20190624'};
selectedFiles = {
    [5:8 14:17]; ...This brings problems
    [1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 20 ];...This  is too much recovered yet
}% each cell is a date
stimMask = {
    [2 2 2 2 2 2 2 2 ]; ...
    [0 0 1 1 0 1 1 0 1 1 1 0 1 1 1 1 0 1]; ... 
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
            VICONfile = [filenameVICON, num2str(iFile, '%02.f')];
            VICONdat{iif} = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '.mat' ]));
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



%%  Preprocessing stage:

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
        % check for complex values
        complexvaluesidx = find(imag(VICONdat{iFile}.jointAngles.jointangles.Elbow) ~= 0);
        VICONdat{iFile}.jointAngles.jointangles.Elbow(complexvaluesidx) = 0;
       
        % Retrieve analogs
        VICONdat{iFile} = getVICONAnalogIn(VICONdat{iFile}, analogs);

        %Retrieve KUKA position
        %VICONdat{iFile} = getKUKAposition(VICONdat{iFile});

        % Syncronization
        % Number of square peaks in TDT trig in VICON
        [VICONdat{iFile}] = alignTDT2VICON_allcases(TDTdat{iFile}, VICONdat{iFile}, 'Brienne', TDTdat{iFile}.date);

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
    if TDTdat{iFile}.isstim == 1 || TDTdat{iFile}.isstim == 0
        [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile},...
            {'startMov', 'grabMov', 'endMov' , 'pull', 'success'});
    else
        [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile},...
            {'startMov', 'grabMov', 'endMov'});
    end
    % Extracting movement attempt indicator
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
    end
    % Assigned success (to do just if you don't filter for success)
    VICONdat{iFile} = assignSuccess(VICONdat{iFile});
    if TDTdat{iFile}.isstim == 2
        VICONdat{iFile} = getTargetPos(VICONdat{iFile});
    end
end



%% Kinematic features computation
        featureslabels = {...
                 'maxSHO', ...% max shoulder angle
                 'maxELB', ...% max elbow angles
                 'maxWRI', ...%max wrist angle
                 'minSHO', ...
                 'minELB', ...
                 'minWRI', ...
                 'excSHO', ...
                 'excELB', ...
                 'excWRI', ...
                 'maxvelSHO', ...
                 'maxvelELB', ...
                 'maxvelWRI', ...
                 'minvelSHO', ...
                 'minvelELB', ...
                 'minvelWRI', ...
                 'avgvelSHO', ...
                 'avgvelELB', ...
                 'avgvelWRI', ...
                 'maxWRIPx', ...
                 'maxWRIPy', ...
                 'maxWRIPz', ...
                 'minWRIPx', ...
                 'minWRIPy', ...
                 'minWRIPz', ...
                 'maxWRIVxreach', ...
                 'maxWRIVyreach', ...
                 'maxWRIVzreach', ...
                 'minWRIVxreach', ...
                 'minWRIVyreach', ...
                 'minWRIVzreach', ...
                 'avgWRIVxreach', ...
                 'avgWRIVyreach', ...
                 'avgWRIVzreach', ...
                 'maxWRIVxpull', ...
                 'maxWRIVypull', ...
                 'maxWRIVzpull', ...
                 'minWRIVxpull', ...
                 'minWRIVypull', ...
                 'minWRIVzpull', ...
                 'avgWRIVxpull', ...
                 'avgWRIVypull', ...
                 'avgWRIVzpull', ...
                 'maxWRIV', ...
                 'minWRIV', ...
                 'avgWRIV', ...
                 'grabWRIPy', ...
                 'grabWRIPz', ...
                 'maxELBPx_R' , ...% >>>>>>>>> TRAJ ELBOW
                 'minELBPx_R' , ...
                 'avgELBPx_R' , ...
                 'maxELBPx_P' , ...
                 'minELBPx_P' , ...
                 'avgELBPx_P' , ...
                 'maxELBPy_R' , ...
                 'minELBPy_R' , ...
                 'avgELBPy_R' , ...
                 'maxELBPy_P' , ...
                 'minELBPy_P' , ...
                 'avgELBPy_P' , ...
                 'maxELBPz_R' , ...
                 'minELBPz_R' , ...
                 'avgELBPz_R' , ...
                 'maxELBPz_P' , ...
                 'minELBPz_P' , ...
                 'avgELBPz_P' , ...
                 'grabELBPx', ...
                 'grabELBPy', ...
                 'grabELBPz', ...
                 'maxELBVx_R' , ...% >>>>>>>>> TRAJ ELBOW
                 'minELBVx_R' , ...
                 'avgELBVx_R' , ...
                 'maxELBVx_P' , ...
                 'minELBVx_P' , ...
                 'avgELBVx_P' , ...
                 'maxELBVy_R' , ...
                 'minELBVy_R' , ...
                 'avgELBVy_R' , ...
                 'maxELBVy_P' , ...
                 'minELBVy_P' , ...
                 'avgELBVy_P' , ...
                 'maxELBVz_R' , ...
                 'minELBVz_R' , ...
                 'avgELBVz_R' , ...
                 'maxELBVz_P' , ...
                 'minELBVz_P' , ...
                 'avgELBVz_P' , ...
                 'maxSHOPx_R' , ...
                 'minSHOPx_R' , ...
                 'avgSHOPx_R' , ...
                 'maxSHOPx_P' , ...
                 'minSHOPx_P' , ...
                 'avgSHOPx_P' , ...
                 'maxSHOPy_R' , ...
                 'minSHOPy_R' , ...
                 'avgSHOPy_R' , ...
                 'maxSHOPy_P' , ...
                 'minSHOPy_P' , ...
                 'avgSHOPy_P' , ...
                 'maxSHOPz_R' , ...
                 'minSHOPz_R' , ...
                 'avgSHOPz_R' , ...
                 'maxSHOPz_P' , ...
                 'minSHOPz_P' , ...
                 'avgSHOPz_P' , ...
                 'grabSHOPy', ...
                 'grabSHOPz', ...
                 'smoothnessR', ...
                 'smoothnessP', ...
                 'maxFx', ...
                 'maxFy', ...
                 'maxFz', ...
                 'maxF', ...
                 'timemaxFx', ...
                 'timeREACH', ...%
                 'timePULL', ...%
                 'lengthREACH', ...%
                 'lengthPULL', ...%
                 };        
fsa = 1000;
fsk = 100;
features = [];
isstim = [];
isstimtimedreach = [];
isstimtimedpull = [];
issuccess = [];
snips = {};
for iFile = 1 : length(VICONdat)
   
    % Filtering
    VICONdat{iFile}.jointAngles.jointangles.Shoulder(isnan(VICONdat{iFile}.jointAngles.jointangles.Shoulder)) = 0;
    VICONdat{iFile}.jointAngles.jointangles.Elbow(isnan(VICONdat{iFile}.jointAngles.jointangles.Elbow)) = 0;
    VICONdat{iFile}.jointAngles.jointangles.Wrist(isnan(VICONdat{iFile}.jointAngles.jointangles.Wrist)) = 0;
    VICONdat{iFile}.jointAngles.wristpos.x(isnan(VICONdat{iFile}.jointAngles.wristpos.x)) = 0;
    VICONdat{iFile}.jointAngles.wristpos.y(isnan(VICONdat{iFile}.jointAngles.wristpos.y)) = 0;
    VICONdat{iFile}.jointAngles.wristpos.z(isnan(VICONdat{iFile}.jointAngles.wristpos.z)) = 0;
    filtSHO = lowPass(VICONdat{iFile}.jointAngles.jointangles.Shoulder', 100, 6, 3);
    filtELB = lowPass(VICONdat{iFile}.jointAngles.jointangles.Elbow', 100, 6, 3);
    filtWRI = lowPass(VICONdat{iFile}.jointAngles.jointangles.Wrist', 100, 6, 3);
    filtvelSHO = diff(filtSHO);
    filtvelELB = diff(filtELB);
    filtvelWRI = diff(filtWRI);
    filtWRIPx = lowPass(VICONdat{iFile}.jointAngles.wristpos.x, 100, 6, 3);
    filtWRIPy= lowPass(VICONdat{iFile}.jointAngles.wristpos.y, 100, 6, 3);
    filtWRIPz= lowPass(VICONdat{iFile}.jointAngles.wristpos.z, 100, 6, 3);
    filtWRIVx = diff(filtWRIPx);
    filtWRIVy = diff(filtWRIPy);
    filtWRIVz = diff(filtWRIPz);
    filtELBx = lowPass(VICONdat{iFile}.kinematic.x(:,3), 100, 6, 3);
    filtELBy = lowPass(-VICONdat{iFile}.kinematic.y(:,3), 100, 6, 3);
    filtELBz = lowPass(VICONdat{iFile}.kinematic.z(:,3), 100, 6, 3);
    
    filtELBVx = diff(filtELBx);
    filtELBVy = diff(filtELBy);
    filtELBVz = diff(filtELBz);
    
    filtSHOx = lowPass(VICONdat{iFile}.kinematic.x(:,1), 100, 6, 3);
    filtSHOy = lowPass(VICONdat{iFile}.kinematic.y(:,1), 100, 6, 3);
    filtSHOz = lowPass(VICONdat{iFile}.kinematic.z(:,1), 100, 6, 3);
    
    filtFx = lowPass(VICONdat{iFile}.force(:,1), 1000, 6, 3);
    filtFy = lowPass(VICONdat{iFile}.force(:,2), 1000, 6, 3);
    filtFz = lowPass(VICONdat{iFile}.force(:,3), 1000, 6, 3);
    if TDTdat{iFile}.isstim == 2
        movcentraltarget = find(VICONdat{iFile}.targetPos == 1);
    else
        movcentraltarget = [1 : length(VICONdat{iFile}.events_tokeep.startMov)];
    end
    for iM =movcentraltarget
        if TDTdat{iFile}.isstim == 2 && not(VICONdat{iFile}.targetPos(iM) ==1)
            continue;
        end
        % Start : present everywhere
        startM = VICONdat{iFile}.events_tokeep.startMov(iM);
        startM_ia = floor(startM*fsa);
        startM_ik = floor(startM*fsk);
        % End : present everywhere 
        endM = VICONdat{iFile}.events_tokeep.endMov(iM);
        endM_ia = floor(endM*fsa);
        endM_ik = floor(endM*fsk);
        % PUll and success : not always presetn
        if TDTdat{iFile}.isstim == 1 || TDTdat{iFile}.isstim == 0
            pullM = VICONdat{iFile}.events_tokeep.pull(iM);
            pullM_ia = floor(pullM*fsa);
            pullM_ik = floor(pullM*fsk);
            grabM = VICONdat{iFile}.events_tokeep.grabMov(iM);
            grabM_ia = floor(grabM*fsa);
            grabM_ik = floor(grabM*fsk);
            mystart_k =startM_ik;  myend_k =endM_ik; mygrab_k = grabM_ik; 
            mystart_a =startM_ia;  myend_a =endM_ia; mygrab_a =grabM_ia;
            timeREACH = grabM - startM;
            timePULL = endM - grabM;
%             successM = VICONdat{iFile}.events_tokeep.success(iM);
%             successM_ia = floor(successM*fsa);
%             successM_ik = floor(successM*fsk);
        else  % grab: not always presetn
            grabM = VICONdat{iFile}.events_tokeep.grabMov(iM);
            grabM_ia = floor(grabM*fsa);
            grabM_ik = floor(grabM*fsk);
            mystart_k =startM_ik;  myend_k =endM_ik; mygrab_k = grabM_ik; 
            mystart_a =startM_ia;  myend_a =endM_ia; mygrab_a = grabM_ia;
            timeREACH = grabM - startM;
            timePULL = endM - grabM;
        end
        
        % Angles
        maxSHO = max(filtSHO(mystart_k:myend_k));
        maxELB = max(filtELB(mystart_k:myend_k));
        maxWRI = max(filtWRI(mystart_k:myend_k));
        minSHO = min(filtSHO(mystart_k:myend_k));
        minELB = min(filtELB(mystart_k:myend_k));
        minWRI = min(filtWRI(mystart_k:myend_k));
        excSHO = maxSHO - minSHO;
        excELB = maxELB - minELB;
        excWRI = maxWRI - minWRI;
        
        maxvelSHO = max(filtSHO(mystart_k:myend_k));
        maxvelELB = max(filtvelELB(mystart_k:myend_k));
        maxvelWRI = max(filtvelWRI(mystart_k:myend_k));
        minvelSHO = max(filtSHO(mystart_k:myend_k));
        minvelELB = max(filtvelELB(mystart_k:myend_k));
        minvelWRI = max(filtvelWRI(mystart_k:myend_k));
        avgvelSHO = mean(filtvelSHO(mystart_k:myend_k));
        avgvelELB = mean(filtvelELB(mystart_k:myend_k));
        avgvelWRI = mean(filtvelWRI(mystart_k:myend_k));
        
        % Wrist trajectory
        maxWRIPx = max(filtWRIPx(mystart_k:myend_k) - filtWRIPx(mystart_k));
        maxWRIPy = max(filtWRIPy(mystart_k:myend_k) - filtWRIPy(mystart_k));
        maxWRIPz = max(filtWRIPz(mystart_k:myend_k) - filtWRIPz(mystart_k));
        minWRIPx = max(filtWRIPx(mystart_k:myend_k) - filtWRIPx(mystart_k));
        minWRIPy = max(filtWRIPy(mystart_k:myend_k) - filtWRIPy(mystart_k));
        minWRIPz = min(filtWRIPz(mystart_k:myend_k) - filtWRIPz(mystart_k));
        
        maxWRIVx_R = max(filtWRIVx(mystart_k:mygrab_k));
        maxWRIVy_R = max(filtWRIVy(mystart_k:mygrab_k));
        maxWRIVz_R = max(filtWRIVz(mystart_k:mygrab_k));
        minWRIVx_R = max(filtWRIVx(mystart_k:mygrab_k));
        minWRIVy_R = max(filtWRIVy(mystart_k:mygrab_k));
        minWRIVz_R = max(filtWRIVz(mystart_k:mygrab_k));
        avgWRIVx_R = mean(filtWRIVx(mystart_k:mygrab_k));
        avgWRIVy_R = mean(filtWRIVy(mystart_k:mygrab_k));
        avgWRIVz_R = mean(filtWRIVz(mystart_k:mygrab_k)); 
        
        maxWRIVx_P = max(filtWRIVx(mygrab_k:myend_k));
        maxWRIVy_P = max(filtWRIVy(mygrab_k:myend_k));
        maxWRIVz_P = max(filtWRIVz(mygrab_k:myend_k));
        minWRIVx_P = max(filtWRIVx(mygrab_k:myend_k));
        minWRIVy_P = max(filtWRIVy(mygrab_k:myend_k));
        minWRIVz_P = max(filtWRIVz(mygrab_k:myend_k));
        avgWRIVx_P = mean(filtWRIVx(mygrab_k:myend_k));
        avgWRIVy_P = mean(filtWRIVy(mygrab_k:myend_k));
        avgWRIVz_P = mean(filtWRIVz(mygrab_k:myend_k)); 
        maxWRIV = max(sqrt(filtWRIVx.^2 + filtWRIVy.^2 + filtWRIVz.^2));
        minWRIV = min(sqrt(filtWRIVx.^2 + filtWRIVy.^2 + filtWRIVz.^2));
        avgWRIV = min(sqrt(filtWRIVx.^2 + filtWRIVy.^2 + filtWRIVz.^2));
        grabWRIPy = mean(filtWRIPy(floor((mygrab_k + myend_k)/2)));
        grabWRIPz = mean(filtWRIPz(floor((mygrab_k + myend_k)/2)));
        
        maxELBPx_R = max(filtELBx(mystart_k:mygrab_k) - filtELBx(mystart_k));
        minELBPx_R = min(filtELBx(mystart_k:mygrab_k)- filtELBx(mystart_k));
        avgELBPx_R = mean(filtELBx(mystart_k:mygrab_k)- filtELBx(mystart_k));
        maxELBPx_P = max(filtELBx(mygrab_k:myend_k) - filtELBx(mygrab_k));
        minELBPx_P = min(filtELBx(mygrab_k:myend_k)- filtELBx(mygrab_k));
        avgELBPx_P = mean(filtELBx(mygrab_k:myend_k)- filtELBx(mygrab_k)); 
        maxELBPy_R = max(filtELBy(mystart_k:mygrab_k) - filtELBy(mystart_k));
        minELBPy_R = min(filtELBy(mystart_k:mygrab_k)- filtELBy(mystart_k));
        avgELBPy_R = mean(filtELBy(mystart_k:mygrab_k)- filtELBy(mystart_k));
        maxELBPy_P = max(filtELBy(mygrab_k:myend_k) - filtELBy(mygrab_k));
        minELBPy_P = min(filtELBy(mygrab_k:myend_k)- filtELBy(mygrab_k));
        avgELBPy_P = mean(filtELBy(mygrab_k:myend_k)- filtELBy(mygrab_k)); 
        maxELBPz_R = max(filtELBz(mystart_k:mygrab_k) - filtELBz(mystart_k));
        minELBPz_R = min(filtELBz(mystart_k:mygrab_k)- filtELBz(mystart_k));
        avgELBPz_R = mean(filtELBz(mystart_k:mygrab_k)- filtELBz(mystart_k));
        maxELBPz_P = max(filtELBz(mygrab_k:myend_k) - filtELBz(mygrab_k));
        minELBPz_P = min(filtELBz(mygrab_k:myend_k)- filtELBz(mygrab_k));
        avgELBPz_P = mean(filtELBz(mygrab_k:myend_k)- filtELBz(mygrab_k));
        grabELBPx = mean(filtELBx(mygrab_k : mygrab_k+ 10)- filtELBx(mystart_k)) ;
        grabELBPy = mean(filtELBy(mygrab_k : mygrab_k+ 10)- filtELBy(mystart_k)) ;
        grabELBPz = mean(filtELBz(mygrab_k : mygrab_k + 10)- filtELBz(mystart_k)) ;
        
        maxELBVx_R = max(filtELBVx(mystart_k:mygrab_k));
        minELBVx_R = min(filtELBVx(mystart_k:mygrab_k));
        avgELBVx_R = mean(filtELBVx(mystart_k:mygrab_k));
        maxELBVx_P = max(filtELBVx(mygrab_k:myend_k));
        minELBVx_P = min(filtELBVx(mygrab_k:myend_k));
        avgELBVx_P = mean(filtELBVx(mygrab_k:myend_k)); 
        maxELBVy_R = max(filtELBVy(mystart_k:mygrab_k));
        minELBVy_R = min(filtELBVy(mystart_k:mygrab_k));
        avgELBVy_R = mean(filtELBVy(mystart_k:mygrab_k));
        maxELBVy_P = max(filtELBVy(mygrab_k:myend_k) );
        minELBVy_P = min(filtELBVy(mygrab_k:myend_k));
        avgELBVy_P = mean(filtELBy(mygrab_k:myend_k)); 
        maxELBVz_R = max(filtELBVz(mystart_k:mygrab_k) );
        minELBVz_R = min(filtELBVz(mystart_k:mygrab_k));
        avgELBVz_R = mean(filtELBVz(mystart_k:mygrab_k));
        maxELBVz_P = max(filtELBVz(mygrab_k:myend_k));
        minELBVz_P = min(filtELBVz(mygrab_k:myend_k));
        avgELBVz_P = mean(filtELBVz(mygrab_k:myend_k));
        
        
        
        
        
        
        maxSHOPx_R = max(filtSHOx(mystart_k:mygrab_k) - filtSHOx(mystart_k));
        minSHOPx_R = min(filtSHOx(mystart_k:mygrab_k) - filtSHOx(mystart_k));
        avgSHOPx_R = mean(filtSHOx(mystart_k:mygrab_k) - filtSHOx(mystart_k));
        maxSHOPx_P = max(filtSHOx(mygrab_k: myend_k) - filtSHOx(mygrab_k));
        minSHOPx_P = min(filtSHOx(mygrab_k:myend_k) - filtSHOx(mygrab_k));
        avgSHOPx_P = mean(filtSHOx(mygrab_k:myend_k) - filtSHOx(mygrab_k));
        maxSHOPy_R = max(filtSHOy(mystart_k:mygrab_k) - filtSHOy(mystart_k));
        minSHOPy_R = min(filtSHOy(mystart_k:mygrab_k) - filtSHOy(mystart_k));
        avgSHOPy_R = mean(filtSHOy(mystart_k:mygrab_k) - filtSHOy(mystart_k));
        maxSHOPy_P = max(filtSHOy(mygrab_k:myend_k) - filtSHOy(mygrab_k));
        minSHOPy_P = min(filtSHOy(mygrab_k:myend_k) - filtSHOy(mygrab_k));
        avgSHOPy_P = mean(filtSHOy(mygrab_k:myend_k) - filtSHOy(mygrab_k));
        maxSHOPz_R = max(filtSHOz(mystart_k:mygrab_k) - filtSHOz(mystart_k));
        minSHOPz_R = min(filtSHOz(mystart_k:mygrab_k) - filtSHOz(mystart_k));
        avgSHOPz_R = mean(filtSHOz(mystart_k:mygrab_k) - filtSHOz(mystart_k));
        maxSHOPz_P = max(filtSHOz(mygrab_k:myend_k) - filtSHOz(mygrab_k));
        minSHOPz_P = min(filtSHOz(mygrab_k:myend_k) - filtSHOz(mygrab_k));
        avgSHOPz_P = mean(filtSHOz(mygrab_k:myend_k) - filtSHOz(mygrab_k));
        grabSHOPy = mean(filtSHOy(floor((mygrab_k + myend_k)/2)));
        grabSHOPz = mean(filtSHOz(floor((mygrab_k + myend_k)/2)));
        
   
        xwrist = lowPass(VICONdat{iFile}.jointAngles.wristpos.x(mystart_k:mygrab_k), 100, 6, 3);
        ywrist = lowPass(VICONdat{iFile}.jointAngles.wristpos.y(mystart_k:mygrab_k), 100, 6, 3);
        zwrist = lowPass(VICONdat{iFile}.jointAngles.wristpos.z(mystart_k:mygrab_k), 100, 6, 3);
        v = sqrt(diff(xwrist).^2 + diff(ywrist).^2 + diff(zwrist).^2);
        time = length(v).*(1/100);
        smoothnessR = length(findpeaks(v))./time ;
        xwrist = lowPass(VICONdat{iFile}.jointAngles.wristpos.x(mygrab_k:myend_k), 100, 6, 3);
        ywrist = lowPass(VICONdat{iFile}.jointAngles.wristpos.y(mygrab_k:myend_k), 100, 6, 3);
        zwrist = lowPass(VICONdat{iFile}.jointAngles.wristpos.z(mygrab_k:myend_k), 100, 6, 3);
        v = sqrt(diff(xwrist).^2 + diff(ywrist).^2 + diff(zwrist).^2);
        time = length(v).*(1/100);
        smoothnessP = length(findpeaks(v))./time;
        %
        
        % Force during pull
        [maxFx, timemaxFx] =(min(filtFx(mygrab_a:myend_a)));
        maxFx = abs(maxFx);
       
        [maxFy, timemaxFy] = max(abs(filtFy(mygrab_a:myend_a)));
        [maxFz, timemaxFz] = min((filtFz(mygrab_a:myend_a)));
        maxFz = abs(maxFz);
        timemaxFx = timemaxFx/fsa;
        aucFx = sum(VICONdat{iFile}.force(mygrab_a:myend_a,1));
        aucFy = sum(VICONdat{iFile}.force(mygrab_a:myend_a,1));
        aucFz = sum(VICONdat{iFile}.force(mygrab_a:myend_a,1));
        if isempty(maxFx) maxFx = 0; end
        if isempty(timemaxFx) timemaxFx = 0; end
        if isempty(aucFx) aucFx = 0; end
        
        F = sqrt(filtFx(mygrab_a:myend_a).^2 + filtFy(mygrab_a:myend_a).^2 + filtFz(mygrab_a:myend_a).^2);
        maxF=  max(F);
        % Movement trajectories length
        xwrist = lowPass(VICONdat{iFile}.jointAngles.wristpos.x(mystart_k:mygrab_k), 100, 6, 3);
        ywrist = lowPass(VICONdat{iFile}.jointAngles.wristpos.y(mystart_k:mygrab_k), 100, 6, 3);
        zwrist = lowPass(VICONdat{iFile}.jointAngles.wristpos.z(mystart_k:mygrab_k), 100, 6, 3);
        lengthREACH = sum(sqrt(diff(xwrist).^2 + diff(ywrist).^2 + diff(zwrist).^2));
        xwrist = VICONdat{iFile}.jointAngles.wristpos.x(mygrab_k:myend_k);
        ywrist = VICONdat{iFile}.jointAngles.wristpos.y(mygrab_k:myend_k);
        zwrist = VICONdat{iFile}.jointAngles.wristpos.z(mygrab_k:myend_k);
        lengthPULL = sum(sqrt(diff(xwrist).^2 + diff(ywrist).^2 + diff(zwrist).^2));
  
       featvect = [maxSHO, ... % >>>>>>>>> ANGLES
         maxELB, ...
         maxWRI, ...
         minSHO, ...
         minELB, ...
         minWRI, ...
         excSHO, ...
         excELB, ...
         excWRI, ...
         maxvelSHO, ...
         maxvelELB, ...
         maxvelWRI, ...
         minvelSHO, ...
         minvelELB, ...
         minvelWRI, ...
         avgvelSHO, ...
         avgvelELB, ...
         avgvelWRI, ... 
         maxWRIPx, ...% >>>>>>>>> TRAJ WRIST
         maxWRIPy, ...
         maxWRIPz, ...
         minWRIPx, ...
         minWRIPy, ...
         minWRIPz, ...
         maxWRIVx_R, ...
         maxWRIVy_R, ...
         maxWRIVz_R, ...
         minWRIVx_R, ...
         minWRIVy_R, ...
         minWRIVz_R, ...
         avgWRIVx_R, ...
         avgWRIVy_R, ...
         avgWRIVz_R, ... 
         maxWRIVx_P, ...
         maxWRIVy_P, ...
         maxWRIVz_P, ...
         minWRIVx_P, ...
         minWRIVy_P, ...
         minWRIVz_P, ...
         avgWRIVx_P, ...
         avgWRIVy_P, ...
         avgWRIVz_P, ...
         maxWRIV, ...
         minWRIV, ...
         avgWRIV, ...
         grabWRIPy, ...
         grabWRIPz, ...
         maxELBPx_R , ...% >>>>>>>>> TRAJ ELBOW
         minELBPx_R , ...
         avgELBPx_R , ...
         maxELBPx_P , ...
         minELBPx_P , ...
         avgELBPx_P , ...
         maxELBPy_R , ...
         minELBPy_R , ...
         avgELBPy_R , ...
         maxELBPy_P , ...
         minELBPy_P , ...
         avgELBPy_P , ...
         maxELBPz_R , ...
         minELBPz_R , ...
         avgELBPz_R , ...
         maxELBPz_P , ...
         minELBPz_P , ...
         avgELBPz_P , ...
         grabELBPx, ...
         grabELBPy, ...
         grabELBPz, ...
         maxELBVx_R , ...% >>>>>>>>> TRAJ VEL ELBOW
         minELBVx_R , ...
         avgELBVx_R , ...
         maxELBVx_P , ...
         minELBVx_P , ...
         avgELBVx_P , ...
         maxELBVy_R , ...
         minELBVy_R , ...
         avgELBVy_R , ...
         maxELBVy_P , ...
         minELBVy_P , ...
         avgELBVy_P , ...
         maxELBVz_R , ...
         minELBVz_R , ...
         avgELBVz_R , ...
         maxELBVz_P , ...
         minELBVz_P , ...
         avgELBVz_P , ...       
         maxSHOPx_R , ...% >>>>>>>>> TRAJ SHOULDER
         minSHOPx_R , ...
         avgSHOPx_R , ...
         maxSHOPx_P , ...
         minSHOPx_P , ...
         avgSHOPx_P , ...
         maxSHOPy_R , ...
         minSHOPy_R , ...
         avgSHOPy_R , ...
         maxSHOPy_P , ...
         minSHOPy_P , ...
         avgSHOPy_P , ...
         maxSHOPz_R , ...
         minSHOPz_R , ...
         avgSHOPz_R , ...
         maxSHOPz_P , ...
         minSHOPz_P , ...
         avgSHOPz_P , ...
         grabSHOPy, ...
         grabSHOPz, ...
         smoothnessR, ...
         smoothnessP, ...
         maxFx, ...
         maxFy, ...
         maxFz, ...
         maxF, ...
         timemaxFx, ...
         timeREACH, ...%
         timePULL, ...%
         lengthREACH, ...%
         lengthPULL];
         if not(isreal(featvect))  
             disp([num2str(iFile),', ', num2str(iM), ' is not real']);
             continue;
         end
        
         if TDTdat{iFile}.isstim == 1
             if TDTdat{iFile}.isstimtimed.pull(iM) == 0 &&  TDTdat{iFile}.isstimtimed.reach(iM) == 0
                 continue;
             end
         end
         features = [features; featvect];
         
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
         
    end  
end

i_eliminate = find(isnan(isstim));
features(i_eliminate, :) = [];
isstim(i_eliminate, :) = [];
isstimtimed(i_eliminate, :) = [];

% barplots of specifi features

nFeat = size(features, 2);
dims = numSubplots(nFeat);
% Compute centroids 
BASidx = find(isstim == 2);
STIMPULLidx = intersect(find(isstim == 1), find(isstimtimedpull == 1));
LESIONidx = find(isstim == 0);




%% Violin plot
iF = 25; % Wrist velocity

figure
dims = numSubplots(length(featureslabels)); 
for  iF = 25%1 : length(featureslabels)
    subplot(dims(1), dims(2), iF)
    baseline = features(BASidx, iF);
    lesion = features(LESIONidx, iF);
    stim = features(STIMPULLidx, iF);

    data = [baseline; lesion; stim];
    datalabels = [ones(size(baseline)); ...
        2*ones(size(lesion)); ...
        3*ones(size(stim))];
    
    %[P,ANOVATAB,STATS] = kruskalwallis(data, datalabels);
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
end


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

subplot(2 ,2, 2)
hold on
scspapercolors = [230 230 230; 125 125 125; 255 183 27];
countbastrials = 0;
for ip = 1 : size(SCORE, 1)
    if issuccess(ip) == 1
        marker = 'v';
        
    else
        marker = 'o';
    end
    if isstim(ip) == 1
        if isstimtimedpull(ip) == 1 || isstimtimedreach(ip) == 1
            if issuccess(ip)
                color = [255 183 27]./255; 
            else
                color = [255 0 27]./255;
            end
        elseif  isstimtimedreach(ip)  == 1
             color = [244 112 22]./255;
             %continue;
        end
    elseif isstim(ip) == 2
        color = [0 0 0];
    else 
        countbastrials = countbastrials + 1;
         color = [125 125 125]./255;
    end
    plot3(SCORE(ip,1), SCORE(ip,2) ,SCORE(ip,3), marker,...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color, ...
        'Markersize',6)

end
% Compute centroids 
BASidx = find(isstim == 2);
STIMREACHidx =  intersect(find(isstim == 1), find(isstimtimedreach == 1));
STIMPULLidx = intersect(find(isstim == 1), find(isstimtimedpull == 1));
LESIONidx = find(isstim == 0);
SUCCESSSTIMidx = intersect(find(isstim == 1), find(issuccess == 1));
SUCCESSBASidx = intersect(find(isstim == 0), find(issuccess == 1));

C_BAS = [mean(SCORE(BASidx, 1)), mean(SCORE(BASidx, 2))] ;
C_STIMPULL = [mean(SCORE(STIMPULLidx, 1)), mean(SCORE(STIMPULLidx, 2))] ;
C_STIMREACH = [mean(SCORE(STIMREACHidx, 1)), mean(SCORE(STIMREACHidx, 2))] ;
C_LESION = [mean(SCORE(LESIONidx, 1)), mean(SCORE(LESIONidx, 2))] ;
C_SUCCESSSTIM = [mean(SCORE(SUCCESSSTIMidx, 1)), mean(SCORE(SUCCESSSTIMidx, 2))] ;
C_SUCCESSBAS = [mean(SCORE(SUCCESSBASidx, 1)), mean(SCORE(SUCCESSBASidx, 2))] ;
hold on
plot(C_BAS(1), C_BAS(2), 'o', 'MarkerFaceColor', [0.3 0.3 0.3],'MarkerEdgeColor', [1 0 0], 'Markersize', 14)
plot(C_STIMPULL(1), C_STIMPULL(2),'o', 'MarkerFaceColor', [255 183 27]./255,'MarkerEdgeColor', [255 183 27]./255, 'Markersize', 12)
%plot(C_STIMREACH(1), C_STIMREACH(2),'o', 'MarkerFaceColor', [ 0  255 0]./255,'MarkerEdgeColor',[ 0  255 0]./255, 'Markersize', 12)
plot(C_LESION(1), C_LESION(2),'o', 'MarkerFaceColor', [125 125 125]./255,'MarkerEdgeColor',[125 125 125]./255, 'Markersize', 10)
plot(C_SUCCESSSTIM(1), C_SUCCESSSTIM(2),'v', 'MarkerFaceColor', [255 183 27]./255,'MarkerEdgeColor',[0 0 0]./255, 'Markersize', 10)
%plot(C_SUCCESSBAS(1), C_SUCCESSBAS(2),'v', 'MarkerFaceColor', [125 125 125]./255,'MarkerEdgeColor',[0 0 0]./255, 'Markersize', 10)

xlabel(['PC1'])
ylabel(['PC2'])
zlabel(['PC3'])


% Euclidean distances
% Compute idx 
clear featdistance_STIMLESION featdistance_BASSTIM featdistance_BASLESION
subplot(2, 2, 4)
hold on
BASidx = find(isstim == 2);
STIMREACHidx =  intersect(find(isstim == 1), find(isstimtimedreach == 1));
STIMPULLidx = intersect(find(isstim == 1), find(isstimtimedpull == 1));
STIMALLidx =  find(isstim == 1);
LESIONidx = find(isstim == 0);
SUCCESSSTIMidx = intersect(find(isstim == 1), find(issuccess == 1));
SUCCESSBASidx = intersect(find(isstim == 0), find(issuccess == 1));
featuresfordistance = zscore(features);
for ibas = 1 : length(BASidx)
    
v1 = featuresfordistance(BASidx(ibas), :);
v2 = mean(featuresfordistance(STIMALLidx, :));
featdistance_BASSTIM(ibas) = sqrt(sum((v1 - v2).^2)); 
    
end

for ibas = 1 : length(BASidx)
    
    v1 = featuresfordistance(BASidx(ibas), :);
    v2 = mean(featuresfordistance(LESIONidx, :));
    featdistance_BASLESION(ibas) = sqrt(sum((v1 - v2).^2)); 
    
end

for iles = 1 : length(LESIONidx)
    
    v1 = featuresfordistance(LESIONidx(iles), :);
    v2 = mean(featuresfordistance(STIMALLidx, :));
    featdistance_STIMLESION(iles) = sqrt(sum((v1 - v2).^2)); 
    
end



randdistr = 0.3*rand(length(featdistance_BASSTIM(:)), 1) - 0.3/2;
plot(ones(size(featdistance_BASSTIM(:)))+randdistr , featdistance_BASSTIM(:), 'o')
plot(1 , mean(featdistance_BASSTIM(:)), 'o', ...
    'MarkerFaceColor', [0 0 0], ...
    'MarkerEdgeColor', [0 0 0],... 
    'Markersize', 10)

randdistr = 0.3*rand(length(featdistance_BASLESION(:)), 1) - 0.3/2;
plot(2*ones(size(featdistance_BASLESION(:))) + randdistr, featdistance_BASLESION(:), 'o')
plot(2 , mean(featdistance_BASLESION(:)), 'o', ...
    'MarkerFaceColor', [0 0 0], ...
    'MarkerEdgeColor', [0 0 0],... 
    'Markersize', 10)

randdistr = 0.3*rand(length(featdistance_STIMLESION(:)), 1) -  0.3/2;
plot(3*ones(size(featdistance_STIMLESION(:))) + randdistr, featdistance_STIMLESION(:), 'o')
plot(3 , mean(featdistance_STIMLESION(:)), 'o', ...
    'MarkerFaceColor', [0 0 0], ...
    'MarkerEdgeColor', [0 0 0],... 
    'Markersize', 10)

p = ranksum(featdistance_BASSTIM(:), featdistance_BASLESION(:)); 
title({'Dist of feats from centroids (features zscored)',  ['G1 vs G2 p = ',num2str(p)]} )
ylabel('Euclidean distance (feature space)')
set(gca, 'xtick', [1 2 3], 'xticklabel', {'PRELESION vs STIM', 'PRELESION vs NOSTIM', 'STIM vs NOSTIM'})
xlim([0.6 3.2])






