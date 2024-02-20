warning('off','all')
clear all
clc
close all
% Root = '/Volumes/Capogrosso/Data/Experiments/';
% currentpath = '/Users/Bea/Documents/PhD/MkAnalysis/Code/scripts/finalscspaper';
Root = 'D:\SERVER\Data\Experiments';
currentpath = 'C:\Users\rnelg1\Desktop\Bea\Repos\unifr_cervicalEES\';
currentpath = 'C:\GIT\unifr_cervicalEES';
Root= 'R:\data_raw\primate\SCS_Monkey';
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
selectedFiles = [26 : 31];
freqs = [20 40 50 60 80 100];

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


%% Force preprocessing stage:
peakforce = [];
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
    
    % Initializing
    TDTdat{iFile}.snipForce = [];
    TDTdat{iFile}.peakforce = [];
    TDTdat{iFile}.featsforce.slopeForce = [];
    
    for iBS = 1 : length(burstsStart)
        
        TDTdat{iFile}.snipForce(:,iBS) = force (burstsStart(iBS): burstsStart(iBS) + avgLenghtBurst);
        [maxVF, imaxVF] = max(force (burstsStart(iBS): burstsStart(iBS) + avgLenghtBurst));
        TDTdat{iFile}.featsforce.peakforce(1, iBS) = maxVF;
        TDTdat{iFile}.featsforce.slopeforce(1, iBS) = (maxVF - force(burstsStart(iBS)))/ (imaxVF*(1/24414.06));    
    end
    peakforce= [peakforce;TDTdat{iFile}.featsforce.peakforce' ];
    TDTdat{iFile}.featsforce.meanpeakforce = median(TDTdat{iFile}.featsforce.peakforce);
    TDTdat{iFile}.featsforce.meanslopeforce = median(TDTdat{iFile}.featsforce.slopeforce);
    
end

% Plot those data
figure
hold on
for iiFile = 1 : length(selectedFiles)
    iFile = selectedFiles(iiFile);
    npoints = length(TDTdat{iFile}.featsforce.peakforce);
    % All values
    plot(iiFile*ones(npoints),  TDTdat{iFile}.featsforce.peakforce, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor',[0.5 0.5 0.5])
    % Mean Value
    plot(iiFile,  TDTdat{iFile}.featsforce.meanpeakforce, 'o','MarkerSize', 10,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
    
end
set(gca, 'xtick', [1 : length(selectedFiles)], 'xticklabel', freqs);
xlabel('Frequency [Hz]')
ylabel('Force [N]')
title('Modulation of Peak Force with stim frequency -- 10 points')



figure
hold on
for iiFile = 1 : length(selectedFiles)
    iFile = selectedFiles(iiFile);
    npoints = length(TDTdat{iFile}.featsforce.slopeforce);
    % All values
    plot(iiFile*ones(npoints),  TDTdat{iFile}.featsforce.slopeforce, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor',[0.5 0.5 0.5])
    % Mean Value
    plot(iiFile,  TDTdat{iFile}.featsforce.meanslopeforce, 'o','MarkerSize', 10,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
    
end
set(gca, 'xtick', [1 : length(selectedFiles)], 'xticklabel', freqs);
xlabel('Frequency [Hz]')
ylabel('Slope [N/s]')
title('Modulation of Slope of force with stim frequency -- 10 points')