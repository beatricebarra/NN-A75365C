warning('off','all')
%clear all
clc
%close all

Root = '/Volumes/SanDisk_Bea/SCS_Monkeys/';
currentpath = '/Users/barrab01/Documents/Repos/unifr_cervicalEES';
%load('/Users/Bea/Documents/Repos/unifr_cervicalEES/Structs/onetrueCM.mat')

Animal = 'Brienne';
expDate = '20190702';
recsystems = {'TDT'};
additionalpath = 'TERMINAL/BrienneTerminal/'
set(0, 'DefaultFigureRenderer', 'painters');

% In this case path is different because we have both behav and termianl in
% the same day
dataPath = fullfile(Root, Animal,expDate, additionalpath);

filenameTDT = [expDate,'_BT-'];
filenameVICON = [expDate,'_',Animal, '_'];
filenameBKR = [expDate,'_',Animal, '_'];

Files = scan_files(dataPath);
nFiles = length(Files);
% Params load
selectedFiles = [26 : 31];
freqs = [20 40 50 60 80 100];
% Params no load
selectedFiles = [12 : 17];
freqs = [20 40 50 60 80 100];

selectedFiles = [12 : 17, 26:31];
freqs = [20 40 50 60 80 100 20 40 50 60 80 100];
onlyfreqs = [20 40 50 60 80 100 ];
for iFile = selectedFiles
    
    disp(['Loading File n ', num2str(iFile) ])
    if any(strcmp(recsystems,'TDT')) || any(strcmp(recsystems,'TDTDyn'))
        TDTfile = [filenameTDT, num2str(iFile, '%01.f')];
        TDTdat{iFile} = load(fullfile(dataPath, [TDTfile, '.mat'] ));
        TDTdat{iFile}.fs = 12207.03;
    end
    if any(strcmp(recsystems,'VICON'))
        VICONfile = [filenameVICON, num2str(iFile, '%02.f')];
        VICONdat{iFile} = load(fullfile(Root, Animal,expDate, 'VICON', [VICONfile, '.mat' ]));
    end
    if any(strcmp(recsystems,'BLACKROCK'))
        BKRfile = [filenameBKR,num2str(iFile, '%03.f')];
        nsFile = fullfile(Root, Animal,expDate, 'BLACKROCK', [BKRfile, '.ns5'] );
        nevfile = fullfile(Root, Animal,expDate, 'BLACKROCK', [BKRfile, '.nev'] );
        BKRdat{iFile}.nsData = openNSx('read', nsFile,'p:double', 's:1');
        BKRdat{iFile}.nevData = openNEV(nevfile,'nomat','nosave');
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


analogs = {'BlackRock', 'TDT_Trigger', 'Kuka_Position'};
for iFile = selectedFiles
    disp(['Processing File n ', num2str(iFile) ])
    if any(strcmp(recsystems,'TDT'))
        % Preprocessing EMGs
        TDTdat{iFile} = processEMGBrienne(TDTdat{iFile}, 'blankart'); 
        
       % figh = plotEMG(TDTdat{iFile}, 'envonraw');
        
        
    end
end


%% Computing energy of signal during stim
dt = (1/TDTdat{iFile}.fs);
for iFile = selectedFiles
    disp(['Processing File n ', num2str(iFile) ])
      
    % Loading Trig and Force
    Trig = TDTdat{iFile}.Trig;
    Trig = lowPass(double(Trig), 24414.06, 5, 3);%signal, fS, fCut, order
    Trig(Trig>=0.0001) = 1;Trig(Trig<0.0001) = 0;
    TDTdat{iFile}.stimBurst = Trig;
    force = abs(lowPass(double(TDTdat{iFile}.AnIO(:,5)), 24414.06, 6, 3) - 2.5).*(10/2.5); % check if conversion factor is correct
    TDTdat{iFile}.force = force;
    % Segmenting Bursts
    [burstsSPeaks, burstsStart] = findpeaks(Trig);
    [burstsEPeaks, burstsEnd] = findpeaks(-Trig);
    avgLenghtBurst = floor(1.2* mean(burstsEnd - burstsStart(1:end-1)));
    TDTdat{iFile}.burstsStart = burstsStart;
    TDTdat{iFile}.burstsEnd = burstsEnd;
    TDTdat{iFile}.avgLenghtBurst = avgLenghtBurst;
    for iBS = 1 : length(TDTdat{iFile}.burstsStart)
        bSTART = floor(TDTdat{iFile}.burstsStart(iBS)/2);
        bSTOP = floor(bSTART + TDTdat{iFile}.avgLenghtBurst/2);

        for iC = 1 : size(TDTdat{iFile}.EMGb,2)
            TDTdat{iFile}.EMGenergy(iBS, iC) = sum(dt.*abs(TDTdat{iFile}.EMGb(bSTART : bSTOP, iC)).^2)./length(TDTdat{iFile}.EMGb(bSTART : bSTOP, iC));
        end
    
    end
    
end

%% Plot of signal energy
% Plot those data
% Choose one muscle
%figure

subplot(2, 1, 1)
muscle = 'BIC'; 
iC = find(cellfun(@(x) strcmp(x, muscle), TDTdat{selectedFiles(1)}.muscles));
iC
hold on
g1 = linspace(1, 0, 15)'; 
greyscale = [g1, g1, g1 ];
energy = [];
labelsenergy = [];
for iiFile = 1 : length(selectedFiles)
    
    iFile = selectedFiles(iiFile);
    if iFile >20
        color = [0.5 0.5 0.5];
        barpos = 1;
    else
        color = [1 1 1];
        barpos = -1;
    end
    npoints = length(TDTdat{iFile}.EMGenergy(:, iC));
%     % All values
% DECOMMENT IF YOU WANT COLOR CODE OF BURST ORDER
%     for iBS = 1 : length(TDTdat{iFile}.EMGenergy(:, iC))
%         plot(iiFile,  TDTdat{iFile}.EMGenergy(iBS, iC), 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor',greyscale(iBS,:))
%     end
    currentPOINTS = TDTdat{iFile}.EMGenergy(:, iC);
     randomxcoord = rand(size(currentPOINTS, 1), 1).*2 - mean(rand(size(currentPOINTS, 1), 1).*2);
    xpos = (randomxcoord + barpos+ freqs(iiFile)*ones(size(currentPOINTS, 1), 1))
    % All values
    plot(xpos,  currentPOINTS, 'o', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',color)
    % Mean Value
    plot(freqs(iiFile)+barpos, median(currentPOINTS), 'o','MarkerSize', 10,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
    b = bar(freqs(iiFile)+barpos, median(currentPOINTS));
    b.FaceColor = [1  1 1];
    er =errorbar(freqs(iiFile)+barpos, median(currentPOINTS),-std(currentPOINTS), std(currentPOINTS));
    er.Color = [1 0 0];
    
    energy= [energy;currentPOINTS ];
    labelsenergy = [labelsenergy; freqs(iiFile)*ones(size(currentPOINTS, 1), 1)];
end
set(gca, 'xtick', onlyfreqs, 'xticklabel', onlyfreqs);
%ylim([0 5*10^-10])

xlabel('Frequency [Hz]')
ylabel('Energy of  EMG ')
title('Modulation of EMG energy with stim frequency -- 10 points')


subplot(2, 1, 2)
muscle = 'FDS'; 
iC = find(cellfun(@(x) strcmp(x, muscle), TDTdat{selectedFiles(1)}.muscles));
iC
hold on
g1 = linspace(1, 0, 15)'; 
greyscale = [g1, g1, g1 ];
energy = [];
labelsenergy = [];
for iiFile = 1 : length(selectedFiles)
    iFile = selectedFiles(iiFile);
    npoints = length(TDTdat{iFile}.EMGenergy(:, iC));
    if iFile >20
        color = [0.5 0.5 0.5];
        barpos = 1;
    else
        color = [1 1 1];
        barpos = -1;
    end
%     % All values
% DECOMMENT IF YOU WANT COLOR CODE OF BURST ORDER
%     for iBS = 1 : length(TDTdat{iFile}.EMGenergy(:, iC))
%         plot(iiFile,  TDTdat{iFile}.EMGenergy(iBS, iC), 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor',greyscale(iBS,:))
%     end
    currentPOINTS = TDTdat{iFile}.EMGenergy(:, iC);
    randomxcoord = rand(size(currentPOINTS, 1), 1).*2 - mean(rand(size(currentPOINTS, 1), 1).*2);
    xpos = (randomxcoord + barpos+ freqs(iiFile)*ones(size(currentPOINTS, 1), 1))
    % All values
    plot(xpos,  currentPOINTS, 'o', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',color)
    % Mean Value
    plot(freqs(iiFile)+barpos, median(currentPOINTS), 'o','MarkerSize', 10,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
    b = bar(freqs(iiFile)+barpos, median(currentPOINTS));
    b.FaceColor = [1  1 1];
    er =errorbar(freqs(iiFile)+barpos, median(currentPOINTS),-std(currentPOINTS), std(currentPOINTS));
    er.Color = [1 0 0];
    
    energy= [energy;currentPOINTS ];
    labelsenergy = [labelsenergy; freqs(iiFile)*ones(size(currentPOINTS, 1), 1)];
end
set(gca, 'xtick', onlyfreqs, 'xticklabel', onlyfreqs);
%ylim([0 1.5*10^-7])
xlabel('Frequency [Hz]')
ylabel('Energy of  EMG ')
title('Modulation of EMG energy with stim frequency -- 10 points')

%% Plot of signal energy
% Plot those data
% Choose one muscle
%figure
figure
subplot(2, 1, 1)
muscle = 'BIC'; 
iC = find(cellfun(@(x) strcmp(x, muscle), TDTdat{selectedFiles(1)}.muscles));
iC
hold on
g1 = linspace(1, 0, 15)'; 
greyscale = [g1, g1, g1 ];
energy = [];
labelsenergy = [];
for iiFile = 1 : length(selectedFiles)
    
    iFile = selectedFiles(iiFile);
    if iFile >20
        color = [0.5 0.5 0.5];
    else
        color = [1 1 1];
    end
    npoints = length(TDTdat{iFile}.EMGenergy(:, iC));
%     % All values
% DECOMMENT IF YOU WANT COLOR CODE OF BURST ORDER
%     for iBS = 1 : length(TDTdat{iFile}.EMGenergy(:, iC))
%         plot(iiFile,  TDTdat{iFile}.EMGenergy(iBS, iC), 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor',greyscale(iBS,:))
%     end
    currentPOINTS = TDTdat{iFile}.EMGenergy(:, iC);
    randomxcoord = rand(size(currentPOINTS, 1), 1).*3 - mean(rand(size(currentPOINTS, 1), 1).*3);
    xpos = (randomxcoord + freqs(iiFile)*ones(size(currentPOINTS, 1), 1))
    % All values
    plot(xpos,  currentPOINTS, 'o', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',color)
    % Mean Value
    plot(freqs(iiFile), median(currentPOINTS), 'o','MarkerSize', 10,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
    
    mymu(iiFile) = median(currentPOINTS);
    mystd(iiFile) = std(currentPOINTS);
    
    energy= [energy;currentPOINTS ];
    labelsenergy = [labelsenergy; freqs(iiFile)*ones(size(currentPOINTS, 1), 1)];
end
set(gca, 'xtick', onlyfreqs, 'xticklabel', onlyfreqs);
%ylim([0 5*10^-10])

xlabel('Frequency [Hz]')
ylabel('Energy of  EMG ')
title('Modulation of EMG energy with stim frequency -- 10 points')


subplot(2, 1, 2)
muscle = 'FDS'; 
iC = find(cellfun(@(x) strcmp(x, muscle), TDTdat{selectedFiles(1)}.muscles));
iC
hold on
g1 = linspace(1, 0, 15)'; 
greyscale = [g1, g1, g1 ];
energy = [];
labelsenergy = [];
for iiFile = 1 : length(selectedFiles)
    iFile = selectedFiles(iiFile);
    npoints = length(TDTdat{iFile}.EMGenergy(:, iC));
    if iFile >20
        color = [0.5 0.5 0.5];
    else
        color = [1 1 1];
    end
%     % All values
% DECOMMENT IF YOU WANT COLOR CODE OF BURST ORDER
%     for iBS = 1 : length(TDTdat{iFile}.EMGenergy(:, iC))
%         plot(iiFile,  TDTdat{iFile}.EMGenergy(iBS, iC), 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor',greyscale(iBS,:))
%     end
    currentPOINTS = TDTdat{iFile}.EMGenergy(:, iC);
    randomxcoord = rand(size(currentPOINTS, 1), 1).*3 - mean(rand(size(currentPOINTS, 1), 1).*3);
    xpos = (randomxcoord + freqs(iiFile)*ones(size(currentPOINTS, 1), 1))
    % All values
    plot(xpos,  currentPOINTS, 'o', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',color)
    % Mean Value
    plot(freqs(iiFile), median(currentPOINTS), 'o','MarkerSize', 10,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
    
    
    energy= [energy;currentPOINTS ];
    labelsenergy = [labelsenergy; freqs(iiFile)*ones(size(currentPOINTS, 1), 1)];
end
set(gca, 'xtick', onlyfreqs, 'xticklabel', onlyfreqs);
%ylim([0 1.5*10^-7])
xlabel('Frequency [Hz]')
ylabel('Energy of  EMG ')
title('Modulation of EMG energy with stim frequency -- 10 points')


%% Stats - ranksum
ncomparison = 6; 
filecouples = [12, 26; 13, 27; 14, 28; 15, 29; 16, 30; 17, 31];
for ifc = 1 : size(filecouples, 1)
    muscles = TDTdat{filecouples(ifc, 1)}.muscles;
    for iM = 1 : length(muscles)
        muscle = TDTdat{filecouples(ifc, 1)}.muscles{iM};
        iC = find(cellfun(@(x) strcmp(x, muscle), TDTdat{selectedFiles(1)}.muscles));
        group1 = TDTdat{filecouples(ifc, 1)}.EMGenergy(:, iC);
        group2 = TDTdat{filecouples(ifc, 2)}.EMGenergy(:, iC);
        testresults.(muscle)(ifc)= ranksum(group1, group2)*ncomparison;
    end
end

%% Stats kruscal wallisalldata
alldata = [];
alldatalab = [];
muscle = 'BIC';
iC = find(cellfun(@(x) strcmp(x, muscle), TDTdat{selectedFiles(1)}.muscles));
        
for ifc = 1 : size(filecouples, 1)
    
    alldata = [alldata; TDTdat{filecouples(ifc, 1)}.EMGenergy(:, iC);]
    alldatalab = [alldatalab;...
        ((ifc-1)*2 + 1)*ones(size(TDTdat{filecouples(ifc, 1)}.EMGenergy(:, iC)))];
    alldata = [alldata; TDTdat{filecouples(ifc, 2)}.EMGenergy(:, iC);]
    alldatalab = [alldatalab; ...
        ((ifc-1)*2 + 2)*ones(size(TDTdat{filecouples(ifc, 2)}.EMGenergy(:, iC)))];



end

[P,ANOVATAB,STATS] = kruskalwallis(alldata, alldatalab);
multcompare(STATS)
