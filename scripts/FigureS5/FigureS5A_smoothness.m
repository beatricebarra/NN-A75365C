%% 

%% Extension
%DLCanalysis

DLCfiles_folder = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/labeled_videos';
tracking_folder = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/trackedfiles';

% Reach set 
selDLCfiles = {'C0294RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0295RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0296RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0297RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0298RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0299RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    };


seltrackedfiles = {'C0294_reachEvents.mat'; ...
    'C0295_reachEvents.mat'; ...
    'C0296_reachEvents.mat'; ...
    'C0297_reachEvents.mat'; ...
    'C0298_reachEvents.mat'; ...
    'C0299_reachEvents.mat'; ...
    };


selDLCfile = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/labeled_videos/C0295RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv';
seltrackedfile ='/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/trackedfiles/C0295_reachEvents.mat';
fs = 60;
freqs = [20 40 50 60 80 100 ];




thisfs = 60;
for iF = 1 : length(selDLCfiles)
    
    selDLCfile = fullfile(DLCfiles_folder, selDLCfiles{iF});
    seltrackedfile =  fullfile(tracking_folder, seltrackedfiles{iF});
    kin = importDLCkinematics(selDLCfile, 'Brienne', '20190702')
    kin = computeDLCangles(kin)
    load(seltrackedfile);

    interpframes = floor(event.frames.*length(kin.shoulder.x)/event.numFrames);
    thisfs = fs.*length(kin.shoulder.x)/event.numFrames;
    kin.starts = [];
    kin.ends = [];
    for iframe = 1 : length(interpframes)
       if strcmp(event.labels{iframe}, 'StartTrial')
            kin.starts = [kin.starts; interpframes(iframe)];
        end
        if strcmp(event.labels{iframe}, 'EndTrial')
            kin.ends = [kin.ends; interpframes(iframe)];
        end
    end
    if iF ==6
        kin.starts(end-1:end) = [];
    end
    cut = min(kin.ends-kin.starts);
    
    
    unfiltered_wrist =[ kin.wrist_c.x, kin.wrist_c.y];
    
  
    for imov = 1 : length(kin.ends)
        
       
        
        % Wrist
        wristsnip = unfiltered_wrist(kin.starts(imov): kin.ends(imov), :);

        VICONdat{iF}.smoothness(imov) = computeSmoothness(wristsnip, 1);
    end
   
end
figure
hold on

for iF = 1 :  length(selDLCfiles)
    bar(iF, median(VICONdat{iF}.smoothness))
    errorbar(iF, median(VICONdat{iF}.smoothness), std(VICONdat{iF}.smoothness))

end


%% FLEXION

%DLCanalysis

DLCfiles_folder = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/labeled_videos';
tracking_folder = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/trackedfiles';

% Reach set 
selDLCfiles = {'C0301RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0302RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0303RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0304RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0305RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    'C0306RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
    };


seltrackedfiles = {'C0301_reachEvents.mat'; ...
    'C0302_reachEvents.mat'; ...
    'C0303_reachEvents.mat'; ...
    'C0304_reachEvents.mat'; ...
    'C0305_reachEvents.mat'; ...
    'C0306_reachEvents.mat'; ...
    };


selDLCfile = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/labeled_videos/C0295RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv';
seltrackedfile ='/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/trackedfiles/C0295_reachEvents.mat';
fs = 60;
freqs = [20 40 50 60 80 100 ];

freqsLab = {'20' '40' '50' '60' '80' '100' };



for iF = 1 : length(selDLCfiles)
    
    selDLCfile = fullfile(DLCfiles_folder, selDLCfiles{iF});
    seltrackedfile =  fullfile(tracking_folder, seltrackedfiles{iF});
    kin = importDLCkinematics(selDLCfile, 'Brienne', '20190702')
    kin = computeDLCangles(kin)
    load(seltrackedfile);

    thisfs = fs.*length(kin.shoulder.x)/event.numFrames;
    interpframes = floor(event.frames.*length(kin.shoulder.x)/event.numFrames);
    
    kin.starts = [];
    kin.ends = [];
    for iframe = 1 : length(interpframes)
       if strcmp(event.labels{iframe}, 'StartTrial')
            kin.starts = [kin.starts; interpframes(iframe)];
        end
        if strcmp(event.labels{iframe}, 'EndTrial')
            kin.ends = [kin.ends; interpframes(iframe)];
        end
    end
    
    unfiltered_wrist =[ kin.wrist_c.x, kin.wrist_c.y];
    
  
    for imov = 1 : length(kin.ends)
        
       
        
        % Wrist
        wristsnip = unfiltered_wrist(kin.starts(imov): kin.ends(imov), :);

        VICONdat{iF}.smoothness(imov) = computeSmoothness(wristsnip, 1);
    end
    
end

figure
hold on

for iF = 1 :  length(selDLCfiles)
    bar(iF, median(VICONdat{iF}.smoothness))
    errorbar(iF, median(VICONdat{iF}.smoothness), std(VICONdat{iF}.smoothness))

end