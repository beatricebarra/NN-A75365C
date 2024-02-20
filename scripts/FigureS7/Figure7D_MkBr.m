warning('off','all')
clear all
clc
close all
Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES/Code/scripts/NatNeuro_finalsubmission';
load('/Users/barrab01/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')

Animal = 'Brienne';
% expDates = { '20190608','20190610', '20190611','20190614', '20190624'}; % 
% %expDates = { '20190815','20190816', '20190820'}; % 
% 
% selectedFiles = { 
%     [1 2 3 9 10]; ... 
%     [2 3 16]; ...
%     [1 8 12]; ...
%     [1 3 13]; ...
%     [1 2 5 8 13 18]; ...
% }% each cell is a date
% stimMask = {
%     [0 0 0 0 0]; ...
%     [0 0 0]; ...
%     [0 0 0]; ...
%     [0 0 0]; ...
%     [0 0 0 0 0 0]; ...
%     
% }
% ---------------------------------------
Animal = 'Brienne';
%expDates = {'20190608', '20190611',  '20190614', '20190624'};
expDates = { '20190529', '20190604', '20190608','20190610', '20190611','20190614', '20190624'};
selectedFiles = {
    [5:8 14:17]; ...
    [1 2 3 ];... 
    [1 2 3  9 10 ]; ... 
    [2 3  16 ]; ...
    [1 8 12]; ...
    [1 3  13]; ...
    [1 2  5  8  13  18  ];...This  is too much recovered yet
}% each cell is a date
stimMask = {
    [2 2 2 2 2 2 2 2 ]; ...
    [0 0 0 ]; ...
    [0 0 0  0 0 ]; ...
    [0 0  0 ]; ...
    [0  0  0]; ...
    [0 0  0]; ...
    [0 0 0  0  0  0 ]; ... 
}% each cell is a date, 1 means stim file, 0 is basel

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

%     if any(strcmp(recsystems,'TDT'))
%         % Preprocessing EMGs
%         if TDTdat{iFile}.isstim == 1
%             TDTdat{iFile} = processEMGBrienne(TDTdat{iFile}, 'blankart', 'filter');%'blankart' 
%         else
%             TDTdat{iFile} = processEMGBrienne(TDTdat{iFile}, 'filter');%'blankart' 
%         end
%         %figh = plotEMG(TDTdat{iFile}, 'envonraw')
%         %figh = plotEMG(TDTdat{iFile}, 'envonrawn')
% 
%     end

    if any(strcmp(recsystems,'VICON'))
        %Computation of joint angles
        VICONdat{iFile}.jointAngles = computeJointAngles(VICONdat{iFile});
        % check for complex values
        %complexvaluesidx = find(imag(VICONdat{iFile}.jointAngles.jointangles.Elbow) ~= 0);
        %VICONdat{iFile}.jointAngles.jointangles.Elbow(complexvaluesidx) = 0;
       
        % Retrieve analogs
        VICONdat{iFile} = getVICONAnalogIn(VICONdat{iFile}, analogs);

        %Retrieve KUKA position
        %VICONdat{iFile} = getKUKAposition(VICONdat{iFile});

        % This puts label in vicon files in which the lbles were deleted
        % for some reason
        if iFile == 20 || iFile == 21 || iFile ++ 22
            VICONdat{iFile}.analog.labels = VICONdat{19}.analog.labels; 
        end

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
            {'startMov', 'grabMov', 'endMov' });
    else
        [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile},...
            {'startMov', 'grabMov', 'endMov'});
    end
%     % Extracting movement attempt indicator
%     % Transforming triggers in blackrock in burst triggers
%     if TDTdat{iFile}.isstim == 1
%         BKRdat{iFile} = extractBursts_BKR(BKRdat{iFile});
%         % Clean BKR bursts
%         BKRdat{iFile} = cleanBursts_BKR(BKRdat{iFile}, Animal, TDTdat{iFile}.date);
%     end
    % Clean force signals 
    VICONdat{iFile} = cleanForceMkBr(VICONdat{iFile});
   
%     %Check stim timing
%     if TDTdat{iFile}.isstim == 1
%         isstimtimed = checkStimTiming(BKRdat{iFile}); 
%         TDTdat{iFile}.isstimtimed = isstimtimed;
%     end
    % Assigned success (to do just if you don't filter for success)
    VICONdat{iFile} = assignSuccess(VICONdat{iFile});
    if TDTdat{iFile}.isstim == 2
        VICONdat{iFile} = getTargetPos(VICONdat{iFile});
    end
end


%% Kinematic features computation
           
fsa = 1000;
fsk = 100;
features = [];
isstim = [];
isstimtimedreach = [];
isstimtimedpull = [];
issuccess = [];
snips = {};
globidx = 0;
alldata = cell(length(expDates),1);
for iFile = 1 : length(VICONdat)
   
    
    
    filtFx = lowPass(VICONdat{iFile}.force(:,1), 1000, 6, 3);
    filtFy = lowPass(VICONdat{iFile}.force(:,2), 1000, 6, 3);
    filtFz = lowPass(VICONdat{iFile}.force(:,3), 1000, 6, 3);

    if TDTdat{iFile}.isstim == 2
        movcentraltarget = find(VICONdat{iFile}.targetPos == 1);
    else
        movcentraltarget = [1 : length(VICONdat{iFile}.events_tokeep.startMov)];
    end
    iFile
    movcentraltarget
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
         grabM = VICONdat{iFile}.events_tokeep.grabMov(iM);
        grabM_ia = floor(grabM*fsa);
        grabM_ik = floor(grabM*fsk);
        mygrab_a = grabM_ia; 
        myend_a = endM_ia; 
        
        F = sqrt(filtFx(mygrab_a:myend_a).^2 + filtFy(mygrab_a:myend_a).^2 + filtFz(mygrab_a:myend_a).^2);
        maxF=  max(F);
        globidx = globidx + 1; 
        perf(globidx, 1) =  maxF;
        thedates{globidx} = TDTdat{iFile}.date;
        
        isstim(globidx) = TDTdat{iFile}.isstim;

        idate = find(cellfun(@(x) strcmp(x, TDTdat{iFile}.date), expDates));
        stimidx = isstim(globidx) + 1;
        alldata{idate} = [alldata{idate};maxF]; 
    end
end

%% Write excel file

save_filename = '/Users/barrab01/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/EDITOR ROUND/SourceData/ExtDataFig7D.xlsx'
Sheetname = 'Mk-Br';     
sel_dates = [1, 4, 5, 6, 7]
writecell(alldata([1 4, 5, 6, 7]),save_filename,'Sheet',Sheetname)
{expDates{sel_dates}}


%%

%% Plot force figure
figure
hold on
vectordata = [];
labeldata = [];
for idate = 1: length(expDates)
    idate
   nostimsem = std(alldata{idate, 1});%./sqrt(length(alldata{idate, 1}));
   randvect = 0.2*(rand(size(alldata{idate, 1}))-0.5);
   b = bar(idate*2 - 0.5 , mean(alldata{idate, 1}))
   err = errorbar(idate*2 - 0.5, mean(alldata{idate, 1}), 0, nostimsem)
   plot((idate*2 - 0.5 + randvect),alldata{idate, 1}, 'ok')
   
  
   vectordata = [vectordata; alldata{idate, 1};]; 
   labeldata = [labeldata; ...
       ((idate-1)*2 + 1)*ones(size(alldata{idate, 1}))]
   
  % linecellarray{((idate-1)*2 + 1)} = alldata{idate, 1};
  
end

figure
%violin(linecellarray([1, 3, 4]))
[P, ANOVATAB, STATS] = kruskalwallis(vectordata, labeldata);
multcompare(STATS)



%% Counting performance
nDates = length(expDates);

nReach = cell(3, 1); nReach{1}= zeros(nDates, 1); nReach{2}= zeros(nDates, 1); nReach{3}= zeros(nDates, 1);
nPull = cell(3, 1); nPull{1}= zeros(nDates, 1); nPull{2}= zeros(nDates, 1); nPull{3}= zeros(nDates, 1);
nSuccess = cell(3, 1); nSuccess{1}= zeros(nDates, 1); nSuccess{2}= zeros(nDates, 1); nSuccess{3}= zeros(nDates, 1);



for iFile = 1 : length(VICONdat)
    iD = find(cellfun(@(x) strcmp(x, TDTdat{iFile}.date), expDates));
    if not(isempty(VICONdat{iFile}.event))
        nR = length(find(cellfun(@(x) strcmp(x, 'Reach'), VICONdat{iFile}.event.labels)));
        nP = length(find(cellfun(@(x) strcmp(x, 'Pull'), VICONdat{iFile}.event.labels)));
        nS = length(find(cellfun(@(x) strcmp(x, 'Success'), VICONdat{iFile}.event.labels)));
        isstim = TDTdat{iFile}.isstim;
        nReach{isstim + 1}(iD) = nReach{isstim + 1}(iD) + nR;
        nPull{isstim + 1}(iD) = nPull{isstim + 1}(iD) + nP;
        nSuccess{isstim + 1}(iD) = nSuccess{isstim + 1}(iD) + nS;
    end
end
%%
figure
subplot(3, 1, 1)
hold on
plot(nReach{1}, '-ok')
plot(nReach{2}, '-oy')
set(gca, 'xtick', [1:nDates], 'xticklabel', expDates)
subplot(3, 1, 2)
hold on
plot(nPull{1}, '-ok')
plot(nPull{2}, '-oy')
set(gca, 'xtick', [1:nDates], 'xticklabel', expDates)
subplot(3, 1, 3)
hold on
plot(nSuccess{1}, '-ok')
plot(nSuccess{2}, '-oy')
set(gca, 'xtick', [1:nDates], 'xticklabel', expDates)



% %% Counting performance
% nDates = length(expDates);
% 
% nReach = zeros(nDates, 1);
% nPull = zeros(nDates, 1);
% nSuccess = zeros(nDates, 1);
% 
% for iFile = 1 : length(VICONdat)
%     iD = find(cellfun(@(x) strcmp(x, TDTdat{iFile}.date), expDates));
%     nR = length(find(cellfun(@(x) strcmp(x, 'Reach'), VICONdat{iFile}.event.labels)));
%     nP = length(find(cellfun(@(x) strcmp(x, 'Pull'), VICONdat{iFile}.event.labels)));
%     nS = length(find(cellfun(@(x) strcmp(x, 'Success'), VICONdat{iFile}.event.labels)));
%     nReach(iD) = nReach(iD) + nR;
%     nPull(iD) = nPull(iD) + nP;
%     nSuccess(iD) = nSuccess(iD) + nS;
% end
% 
% figure
% subplot(2, 1, 1)
% hold on
% plot(nReach, '-ob')
% plot(nPull, '-oy')
% plot(nSuccess, '-ok')
% set(gca, 'xtick', [1:nDates], 'xticklabel', expDates)
% 
% subplot(2, 1, 2)
% hold on
% bar([nReach, nPull, nSuccess])
% set(gca, 'xtick', [1:nDates], 'xticklabel', expDates)
