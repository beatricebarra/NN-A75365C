warning('off','all')
clear all
clc
close all
Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES/Code/scripts/NatNeuro_finalsubmission';
load('/Users/barrab01/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')

Animal = 'Brienne';
%expDates = {'20190608', '20190611',  '20190614', '20190624'};
expDates = { '20190529', '20190604', '20190608','20190610', '20190611','20190614', '20190624'};
selectedFiles = {
    [5:8 14:17]; ...
    [1 2 3 ];... 
    [1 2 3  9 10 ]; ... 
    [2 3 16 ]; ...
    [1  8  12]; ...
    [1 3 13]; ...
    [1 2  5  8  13  18 ];...This  is too much recovered yet
}% each cell is a date
stimMask = {
    [2 2 2 2 2 2 2 2 ]; ...
    [0 0 0]; ...
    [0 0 0  0 0 ]; ...
    [0 0 0]; ...
    [0 0  0]; ...
    [0 0 0]; ...
    [0 0  0  0  0  0 ]; ... 
}%
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
        if strcmp(expDate, '20190529') || strcmp(expDate, '20190604')
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

reaches_permin_B = cell(length(expDates), 1); 
reachnum = cell(length(expDates), 1); 
graspnum = cell(length(expDates), 1); 
pullnum = cell(length(expDates), 1); 
totaltime= cell(length(expDates), 1); 
for i = 1 : length(reaches_permin_B)
    reaches_permin_B{i} = [];
    reachnum{i} = 0;
    graspnum{i} = 0;
    pullnum{i} = 0;
    totaltime{i} = 0;
end
for iFile = 1 : length(TDTdat)
    
    iD = find(cellfun(@(x) strcmp(x, TDTdat{iFile}.date), expDates)); 
    
    [VICONdat{iFile}] = alignTDT2VICON_allcases(TDTdat{iFile}, VICONdat{iFile}, 'Brienne', TDTdat{iFile}.date);

     VICONdat{iFile} = getVICONEvents(VICONdat{iFile});
    % Cleaning events

    [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile},...
        {'startMov', 'grabMov', 'endMov'});

        time_trial = size( VICONdat{iFile}.analog.data, 1)*(1/1000);

    isstim(iFile) = TDTdat{iFile}.isstim;  
    trialtimes(iFile) = time_trial;
    %%resample per minute
    nsecs = 60;
    minnum = round(time_trial/nsecs);
    
    good = find(cellfun(@(x) strcmp(x, 'EndTrial'), VICONdat{iFile}.event.labels));
    reaches = find(cellfun(@(x) strcmp(x, 'Reach'), VICONdat{iFile}.event.labels));
    grasps = find(cellfun(@(x) strcmp(x, 'Pull'), VICONdat{iFile}.event.labels));
    pulls = find(cellfun(@(x) strcmp(x, 'Success'), VICONdat{iFile}.event.labels));
    
    if iD == 1
        reachtimes = VICONdat{iFile}.event.times(good); 
        grasptimes =  VICONdat{iFile}.event.times(good); 
        pulltimes = VICONdat{iFile}.event.times(good); 
    else
        reachtimes = VICONdat{iFile}.event.times(reaches);
        grasptimes =  VICONdat{iFile}.event.times(grasps); 
        pulltimes = VICONdat{iFile}.event.times(pulls); 
    end
    reachnum{iD} = reachnum{iD} + length(reachtimes); 
    graspnum{iD} = graspnum{iD} + length(grasptimes); 
    pullnum{iD} = pullnum{iD} + length(pulltimes); 
    totaltime{iD} = totaltime{iD} + time_trial;
    
   
    for imin = 1 : minnum
        if not(imin == minnum)
            reaches_oneminute = intersect(find(reachtimes > (imin-1)*nsecs ), find(reachtimes < (imin)*nsecs ));
            reaches_permin_B{iD} = [reaches_permin_B{iD}, length(reaches_oneminute)]; 

        else
            reaches_oneminute = intersect(find(reachtimes > (imin-1)*nsecs ), find(reachtimes < time_trial ));
            reaches_permin_B{iD} = [reaches_permin_B{iD}, length(reaches_oneminute)]; 


        end
    end
    
end




figure
hold on
subplot(1, 3, 1)
hold on
for iD = 1:length(reachnum)
    b = bar(iD, reachnum{iD}/totaltime{iD}*60);
end 
set(gca, 'xtick', [1:length(reachnum)], 'xticklabel', expDates)
subplot(1, 3, 2)
hold on
for iD = 1:length(reachnum)
    b = bar(iD, graspnum{iD}/totaltime{iD}*60);
end 
set(gca, 'xtick', [1:length(reachnum)], 'xticklabel', expDates)
subplot(1, 3, 3)
hold on
for iD = 1:length(reachnum)
    b = bar(iD, pullnum{iD}/totaltime{iD}*60);
end 
set(gca, 'xtick', [1:length(reachnum)], 'xticklabel', expDates)% figure
figure
hold on
for iD = 1:length(reachnum)
    b = bar(iD, reachnum{iD});
end 
%% 
figure
hold on
for iD = 1:length(reaches_permin_B)
    b = bar(iD, mean(reaches_permin_B{iD}));
    er = errorbar(iD, mean(reaches_permin_B{iD}), 0, std(reaches_permin_B{iD}));
    plot(iD*ones(size(reaches_permin_B{iD})), reaches_permin_B{iD}, 'ok')
end
vect = []; lab = [];
for iD = 1:length(reaches_permin_B)
    vect = [vect ; reaches_permin_B{iD}'];
    lab = [lab ; iD*ones(size(reaches_permin_B{iD}))'];
end
[P, anovatab, stats] = kruskalwallis(vect, lab);
multcompare(stats);