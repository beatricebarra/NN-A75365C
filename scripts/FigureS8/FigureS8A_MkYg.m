warning('off','all')
clear all
clc
close all
Root = 'P:\data_raw\primate\SCS_Monkey\';
currentpath = 'C:\GIT\unifr_cervicalEES\';
load('C:\GIT\unifr_cervicalEES\Structs\onetrueCM.mat')


Animal = 'Ygritte';

% expDates = { '20190809', '20190812', '20190813', '20190820'}; % 
% selectedFiles = {
%     [ 4 5 6 7 10 14 15 16 ];...
%     [ 2 8 10  12 13  17 18]; ...% add the 14
%     [2 3 4 10 11 12 13 14 15 16]; ...%2 
%     [1 3 4 6 8 9 10 11 12 14 17]; ...%2    
%     }
% stimMask = {
%     [0 0 1 1 1 1 1 1]; ...
%     [1 1 1 1  1  1 1]; ...
%     [0 0 1 1 1 1 1 1 1 1]; ...%0
%     [0 0 1 1 1 1 1 1 1 1 1]; ...%0  
% }% each cell is a date, 1 means stim file, 0 is baseline


expDates = { '20190813', '20190815' };%

selectedFiles = {
    
   %[ 4 5 10 11 12 14 15 16 ];...% 09-08
   %[ 2  5 13  17 18]; ...% 12-08
   [ 2 3  8 9 11 12 13 14 15 16]; ...%13-08
   [1 2 3 4 5 6 7 8] ; ...% 15-08
   %[1 3 4 6 8 9 10 11 12 14 17]; ...%20 -08
    }
stimMask = {
    
    %[0 0 1 1 1 1 1 1]; ...
    %[ 0 1 1 1 1]
    [ 0 0 2 2 1 1 1 1 1 1]; ...%0
    [0 0 1 1 0 1 1 1 ]; 
    %[0 0 1 1 1 1 1 1 1 1 1]; ...%0  

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


    if any(strcmp(recsystems,'VICON'))


        %Retrieve KUKA position
        %VICONdat{iFile} = getKUKAposition(VICONdat{iFile});

        % Syncronization
        % Number of square peaks in TDT trig in VICON
        %VICONdat{iFile}.TDT_Trigger.times = alignTDT2VICON_robust(TDTdat{iFile}, VICONdat{iFile}, 'Brienne');
        [VICONdat{iFile}] = alignTDT2VICON_allcases(TDTdat{iFile}, VICONdat{iFile}, 'Ygritte', TDTdat{iFile}.date);

        [peak, VICONdat{iFile}.BKR_trigger.samples] = (findpeaks(VICONdat{iFile}.analog.data(:,1), 'MinPeakHeight', 1, 'MinPeakProminence', 1));
        VICONdat{iFile}.BKR_trigger.times = VICONdat{iFile}.BKR_trigger.samples(1)/VICONdat{iFile}.analog.framerate;

        
        
        
    end
    
    if any(strcmp(recsystems,'BLACKROCK'))
        BKRdat{iFile}.TDT_Trigger.times = VICONdat{iFile}.TDT_Trigger.times - VICONdat{iFile}.BKR_trigger.times;
        BKRdat{iFile}.TDT_Trigger.samples = floor(BKRdat{iFile}.TDT_Trigger.times * BKRdat{iFile}.nsData.MetaTags.SamplingFreq);
    end
end
%% Additional processing (more advanced and optional)
% % allmov = cell(2, 1);
% % reachmov= cell(2, 1);
% % reachandpullmov =cell(2, 1);
% % successmov =cell(2, 1);

% for i = 1 : 2
%     allmov{i}= 0;
%     reachmov{i}= 0;
%     reachandpullmov{i} =0;
%     successmov{i} =0;
% end
reaches_permin_B = [];
grasps_permin_B = [];
pulls_permin_B = [];
reaches_permin_S = [];
grasps_permin_S = [];
pulls_permin_S = [];
reaches_permin_T = [];
grasps_permin_T = [];
pulls_permin_T = [];
totalminutesS = 0;
totalminutesB = 0;
totalminutesT = 0;
for iFile = 1 : length(TDTdat)
    
    
    

    time_trial = size( VICONdat{iFile}.analog.data, 1)*(1/1000);

    isstim(iFile) = TDTdat{iFile}.isstim;  
    trialtimes(iFile) = time_trial;
    %%resample per minute
    nsecs = 2;
    minnum = round(time_trial/nsecs);
    
    reaches = find(cellfun(@(x) strcmp(x, 'Reach'), VICONdat{iFile}.event.labels));
    grasps = find(cellfun(@(x) strcmp(x, 'Pull'), VICONdat{iFile}.event.labels));
    pulls = find(cellfun(@(x) strcmp(x, 'Success'), VICONdat{iFile}.event.labels));
    
    
    reachtimes = VICONdat{iFile}.event.times(reaches); 
    grasptimes = VICONdat{iFile}.event.times(grasps); 
    pulltimes = VICONdat{iFile}.event.times(pulls); 
    if  TDTdat{iFile}.isstim == 0
        totalminutesB = totalminutesB + minnum; 
        for imin = 1 : minnum
            if not(imin == minnum)
                reaches_oneminute = intersect(find(reachtimes > (imin-1)*nsecs ), find(reachtimes < (imin)*nsecs ));
                reaches_permin_B = [reaches_permin_B, length(reaches_oneminute)]; 
               
                grasps_oneminute = intersect(find(grasptimes > (imin-1)*nsecs ), find(grasptimes < (imin)*nsecs ));
                grasps_permin_B = [grasps_permin_B, length(grasps_oneminute)]; 
                
                pulls_oneminute = intersect(find(pulltimes > (imin-1)*nsecs ), find(pulltimes < (imin)*nsecs ));
                pulls_permin_B = [pulls_permin_B, length(pulls_oneminute)]; 
            else
                reaches_oneminute = intersect(find(reachtimes > (imin-1)*nsecs ), find(reachtimes < time_trial ));
                reaches_permin_B = [reaches_permin_B, length(reaches_oneminute)]; 
                
                grasps_oneminute = intersect(find(grasptimes > (imin-1)*nsecs ), find(grasptimes < time_trial));
                grasps_permin_B = [grasps_permin_B, length(grasps_oneminute)]; 
                
                pulls_oneminute = intersect(find(pulltimes > (imin-1)*nsecs ), find(pulltimes < time_trial));
                pulls_permin_B = [pulls_permin_B, length(pulls_oneminute)]; 
            end
        end
    elseif  TDTdat{iFile}.isstim == 1
        totalminutesS = totalminutesS + minnum; 
        for imin = 1 : minnum
            if not(imin == minnum)
                reaches_oneminute = intersect(find(reachtimes > (imin-1)*nsecs ), find(reachtimes < (imin)*nsecs ));
                reaches_permin_S = [reaches_permin_S, length(reaches_oneminute)]; 
               
                grasps_oneminute = intersect(find(grasptimes > (imin-1)*nsecs ), find(grasptimes < (imin)*nsecs ));
                grasps_permin_S = [grasps_permin_S, length(grasps_oneminute)]; 
                
                pulls_oneminute = intersect(find(pulltimes > (imin-1)*nsecs ), find(pulltimes < (imin)*nsecs ));
                pulls_permin_S = [pulls_permin_S, length(pulls_oneminute)]; 
            else
                reaches_oneminute = intersect(find(reachtimes > (imin-1)*nsecs ), find(reachtimes < time_trial ));
                reaches_permin_S = [reaches_permin_S, length(reaches_oneminute)]; 
                
                grasps_oneminute = intersect(find(grasptimes > (imin-1)*nsecs ), find(grasptimes < time_trial));
                grasps_permin_S = [grasps_permin_S, length(grasps_oneminute)]; 
                
                pulls_oneminute = intersect(find(pulltimes > (imin-1)*nsecs ), find(pulltimes < time_trial));
                pulls_permin_S = [pulls_permin_S, length(pulls_oneminute)]; 
            end
        end
    elseif TDTdat{iFile}.isstim == 2
        totalminutesT = totalminutesT + minnum; 
        for imin = 1 : minnum
            if not(imin == minnum)
                reaches_oneminute = intersect(find(reachtimes > (imin-1)*nsecs ), find(reachtimes < (imin)*nsecs ));
                reaches_permin_T = [reaches_permin_T, length(reaches_oneminute)]; 
               
                grasps_oneminute = intersect(find(grasptimes > (imin-1)*nsecs ), find(grasptimes < (imin)*nsecs ));
                grasps_permin_T = [grasps_permin_T, length(grasps_oneminute)]; 
                
                pulls_oneminute = intersect(find(pulltimes > (imin-1)*nsecs ), find(pulltimes < (imin)*nsecs ));
                pulls_permin_T = [pulls_permin_T, length(pulls_oneminute)]; 
            else
                reaches_oneminute = intersect(find(reachtimes > (imin-1)*nsecs ), find(reachtimes < time_trial ));
                reaches_permin_T = [reaches_permin_T, length(reaches_oneminute)]; 
                
                grasps_oneminute = intersect(find(grasptimes > (imin-1)*nsecs ), find(grasptimes < time_trial));
                grasps_permin_T = [grasps_permin_T, length(grasps_oneminute)]; 
                
                pulls_oneminute = intersect(find(pulltimes > (imin-1)*nsecs ), find(pulltimes < time_trial));
                pulls_permin_T = [pulls_permin_T, length(pulls_oneminute)]; 
            end
        end
    end
    
    
end


%%
% Reach
figure
N_samples = 10000;
[reaches_boot ]=bootstrap_bea(reaches_permin_B,N_samples);
reaches_dist_B = mean(reaches_boot); 
[reaches_boot ]=bootstrap_bea(reaches_permin_S,N_samples);
reaches_dist_S = mean(reaches_boot); 
[reaches_boot ]=bootstrap_bea(reaches_permin_T,N_samples);
reaches_dist_T = mean(reaches_boot); 
mu = mean(reaches_dist_S); 
reaches_dist_B = reaches_dist_B./mu;
reaches_dist_S = reaches_dist_S./mu;
reaches_dist_T = reaches_dist_T./mu;

M = [mean(reaches_dist_B), mean(reaches_dist_S), mean(reaches_dist_T)]; 
S = [std(reaches_dist_B), std(reaches_dist_S), std(reaches_dist_T) ];
i1 = 2; 
i2 = 3; 
p=cdf('Normal',0,M(i1)-M(i2),sqrt(S(1)^2+S(2)^2))+ cdf('Normal',2*(-M(i1)+M(i2)),-M(i1)+M(i2),sqrt(S(i1)^2+S(i2)^2));

subplot(1, 3, 1)
hold on
b = bar([1, 2, 3],  (M));
er = errorbar([1, 2, 3],  (M), [0 0 0], (S));
title(p)

% Grasp 

[grasps_boot ]=bootstrap_bea(grasps_permin_B,N_samples);
grasps_dist_B = mean(grasps_boot); 
[grasps_boot ]=bootstrap_bea(grasps_permin_S,N_samples);
grasps_dist_S = mean(grasps_boot);
[grasps_boot ]=bootstrap_bea(grasps_permin_T,N_samples);
grasps_dist_T = mean(grasps_boot);

mu = mean(grasps_dist_S); 
grasps_dist_B = grasps_dist_B./mu;
grasps_dist_S = grasps_dist_S./mu;
grasps_dist_T = grasps_dist_T./mu;
M = [mean(grasps_dist_B), mean(grasps_dist_S), mean(grasps_dist_T)]; 
S = [std(grasps_dist_B), std(grasps_dist_S), std(grasps_dist_T)];
%M = fliplr(M);
%S = fliplr(S);
i1 = 2; 
i2 = 1; 
p=cdf('Normal',0,M(i1)-M(i2),sqrt(S(1)^2+S(2)^2))+ cdf('Normal',2*(-M(i1)+M(i2)),-M(i1)+M(i2),sqrt(S(i1)^2+S(i2)^2));
subplot(1, 3, 2)
hold on
b = bar([1, 2, 3],  (M));
er = errorbar([1, 2, 3],  (M), [0 0 0 ], (S));
title(p)
% Pulls 

[pulls_boot ]=bootstrap_bea(pulls_permin_B,N_samples);
pulls_dist_B = mean(pulls_boot); 
[pulls_boot ]=bootstrap_bea(pulls_permin_S,N_samples);
pulls_dist_S = mean(pulls_boot); 
[pulls_boot ]=bootstrap_bea(pulls_permin_T,N_samples);
pulls_dist_T = mean(pulls_boot); 
mu = mean(pulls_dist_S); 
pulls_dist_B = pulls_dist_B./mu;
pulls_dist_S = pulls_dist_S./mu;
pulls_dist_T = pulls_dist_T./mu;
M = [mean(pulls_dist_B), mean(pulls_dist_S), mean(pulls_dist_T)]; 
S = [std(pulls_dist_B), std(pulls_dist_S),  std(pulls_dist_T)];

i1 = 2; 
i2 = 1; 
p=cdf('Normal',0,M(i1)-M(i2),sqrt(S(1)^2+S(2)^2))+ cdf('Normal',2*(-M(i1)+M(i2)),-M(i1)+M(i2),sqrt(S(i1)^2+S(i2)^2));

p 

subplot(1, 3, 3)
hold on
b = bar([1, 2, 3],  (M));
er = errorbar([1, 2, 3],  (M), [0 0 0], (S));
title(p)


