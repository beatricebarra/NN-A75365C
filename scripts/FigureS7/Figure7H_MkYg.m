warning('off','all')
clear all
clc
close all
Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES/Code/scripts/NatNeuro_finalsubmission';
load('/Users/barrab01/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')

Animal = 'Ygritte';
% expDates = { '20190809', '20190812','20190813', '20190820'}; % 
% %expDates = { '20190815','20190816', '20190820'}; % 
% % '20190815','20190816',
% selectedFiles = { 
%     [4 5 ];...
%     [2]; ...
%     [2 3]; ...%2  
%     %[ 2 5]; ...
%     %[1 2 ]; ...
%     [1 3 ]; ...%2    
% }% each cell is a date
% stimMask = {
%     [0 0 ]; ...
%     [0]; ...
%     [0 0 ]; ...%0 
%     %[ 0 0]; ...
%     %[0 0 ]; ...
%     [0 0 ]; ...%0 
% }
% ---------------------------------------

expDates = { '20190802', '20190809', '20190813', '20190820'}; % 
selectedFiles = {
    [ 1 2 ]; ...
    [ 4 5 6 7 10 14 15 16 ];...
   % [ 2 8 10  12 13  17 18]; ...% add the 14
    [2 3 4 10 11 12 13 14 15 16]; ...%2 
    [1 3 4 6 8 9 10 11 12 14 17]; ...%2    
    }
stimMask = {
    [2 2]; ...
    [0 0 1 1 1 1 1 1]; ...
    %[1 1 1 1  1  1 1]; ...
    [0 0 1 1 1 1 1 1 1 1]; ...%0
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
        %VICONdat{iFile}.jointAngles = computeJointAngles(VICONdat{iFile});
        % check for complex values
        %complexvaluesidx = find(imag(VICONdat{iFile}.jointAngles.jointangles.Elbow) ~= 0);
        %VICONdat{iFile}.jointAngles.jointangles.Elbow(complexvaluesidx) = 0;
       
        % Retrieve analogs
        VICONdat{iFile} = getVICONAnalogIn(VICONdat{iFile}, analogs);

        %Retrieve KUKA position
        %VICONdat{iFile} = getKUKAposition(VICONdat{iFile});

        % Syncronization
        % Number of square peaks in TDT trig in VICON
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
    if TDTdat{iFile}.isstim == 1 || TDTdat{iFile}.isstim == 0
        [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile},...
            {'startMov', 'grabMov', 'endMov' });
    else
        [TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile}] = cleanEvents(TDTdat{iFile}, VICONdat{iFile}, BKRdat{iFile},...
            {'startMov', 'grabMov', 'endMov'});
    end
    VICONdat{iFile} = cleanForceMkBr(VICONdat{iFile});
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

    
    movcentraltarget = [1 : length(VICONdat{iFile}.events_tokeep.startMov)];
    
    for iM =movcentraltarget
        
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

%%
save_filename = '/Users/barrab01/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/EDITOR ROUND/SourceData/ExtDataFig7D.xlsx'
Sheetname = 'Mk-Yg';     
sel_dates = [1, 2, 3, 4]
writecell(alldata(sel_dates),save_filename,'Sheet',Sheetname)
{expDates{sel_dates}}