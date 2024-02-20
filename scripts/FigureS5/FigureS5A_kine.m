%DLCanalysis

DLCfiles_folder = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/labeled_videos';
tracking_folder = '/Volumes/DATABACKUP1/SERVER/Data/Experiments/Brienne/20190702/TERMINAL/trackedfiles';

% Extensions set 
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


% % FLexion set 
% selDLCfiles = {'C0301RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
%     'C0302RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
%     'C0303RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
%     'C0304RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
%     'C0305RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
%     'C0306RRCDLC_resnet50_An_MkBrMay7shuffle1_1030000.csv'; ...
%     };
% 
% 
% seltrackedfiles = {'C0301_reachEvents.mat'; ...
%     'C0302_reachEvents.mat'; ...
%     'C0303_reachEvents.mat'; ...
%     'C0304_reachEvents.mat'; ...
%     'C0305_reachEvents.mat'; ...
%     'C0306_reachEvents.mat'; ...
%     };

fs = 60;
freqs = [20 40 50 60 80 100 ];



freqsLab = {'20''40' '60' '50' '80' '100' };



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
    
    if iF == 6
         kin.starts(8) = [];
        kin.starts(end) = [];
    end
        
    cut = min(kin.ends-kin.starts);
    
    filtered_wrist = lowPass(kin.angles.wrist_flex', thisfs, 6, 3);
    filtered_pinkie = lowPass(kin.angles.pinkie_flex', thisfs, 6, 3);
    figure
    for imov = 1 : length(kin.ends)
        
        movlength = kin.ends(imov) - kin.starts(imov);
        
        % Elbow
        elbowsnip = kin.angles.elbow(kin.starts(imov): kin.ends(imov));
        elbowSnips{iF}(:,imov) = elbowsnip(1:cut) - elbowsnip(1);
        [minpeak, minpeakloc]= min(elbowsnip);
        elbowSlope{iF}(imov) = (kin.angles.elbow(kin.starts(imov))- minpeak)/(minpeakloc);
        elbowExcursion{iF}(imov) = elbowsnip(1)- elbowsnip(minpeakloc);
        
        elbowRiseTime{iF}(imov) =(minpeakloc);
        
        elbowMaxDer{iF}(imov) =max(diff(elbowsnip));
        
        % Wrist
        wristsnip = filtered_wrist(kin.starts(imov): kin.ends(imov));
        [maxpeak, maxpeakloc]= max(wristsnip(1 : floor(movlength/2)));
        
%         plot(wristsnip)
%         hold on, plot(maxpeakloc, maxpeak, 'o')
%         pause
        
        wristSnips{iF}(:,imov) = wristsnip(1:cut) - wristsnip(1);
        wristSlope{iF}(imov) = (maxpeak - filtered_wrist(kin.starts(imov)))/(maxpeakloc);
        wristExcursion{iF}(imov) = wristsnip(maxpeakloc)- wristsnip(1);
        wristRiseTime{iF}(imov) =(maxpeakloc);
        wristMaxDer{iF}(imov) =max(diff(wristsnip));
        % Pinkie
        pinkiesnip = filtered_pinkie(kin.starts(imov): kin.ends(imov));
        %[maxpeak, maxpeakloc]= max(pinkiesnip);
        
        [minpeak, minpeakloc]= min(pinkiesnip(1 : floor(movlength/2)));
%         plot(pinkiesnip)
%         hold on, plot(minpeakloc, minpeak, 'o')
%         pause
        pinkieSnips{iF}(:,imov) = pinkiesnip(1:cut) - pinkiesnip(1);
        pinkieSlope{iF}(imov) = (minpeak - filtered_wrist(kin.starts(imov)))/(minpeakloc);
        pinkieExcursion{iF}(imov) = pinkiesnip(minpeakloc)- pinkiesnip(1);
        pinkieRiseTime{iF}(imov) =(minpeakloc);
        pinkieMaxDer{iF}(imov) =max(diff(pinkiesnip));
    end
%     figure
%     plot(kin.angles.pinkie_flex)
%     hold on
%     miny = min(kin.angles.pinkie_flex);
%     maxy = max(kin.angles.pinkie_flex);
% 
%     for i = 1 : length(event.frames)
%         if strcmp(event.labels{i}, 'StartTrial')
%             plot(ones(1, 10).*interpframes(i), linspace(miny, maxy, 10), 'r')
%         end
%         if strcmp(event.labels{i}, 'EndTrial')
%             plot(ones(1, 10).*interpframes(i), linspace(miny, maxy, 10), 'b')
%         end
%     end
end


%% Figure with all snips

figure
jet= colormap(jet);
L = size(jet, 1);
ncols = 6;
colorstep = floor(L/6);


subplot(3, 1, 1)
hold on
for iF = 1 : length(wristExcursion)
    
    stdshade(elbowSnips{iF}(1:30, :)',0.2,jet(colorstep*(iF-1) +1, :))
    
    
    
end
legend(freqsLab)
title('Elbow entension')


subplot(3, 1, 2)
hold on
for iF = 1 : length(wristExcursion)
    
    stdshade(wristSnips{iF}(1:30, :)',0.2,jet(colorstep*(iF-1) +1, :))
    
    
    
end
legend(freqsLab)
title('Wrist flexion')

subplot(3, 1, 3)
hold on
for iF = 1 : length(wristExcursion)
    
    stdshade(pinkieSnips{iF}(1:30, :)',0.2,jet(colorstep*(iF-1) +1, :))
    
    
    
end
legend(freqsLab)
title('Pinkie flexion')


%% Figure on Elbow excursion

figure

subplot(1, 3, 1)
hold on
for iF = 1 : length(elbowExcursion)
    plot(iF*ones(size(elbowExcursion{iF}))+ 0.2*(rand(length(wristSlope{iF}),1)-0.5), elbowExcursion{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(elbowExcursion{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

ylim([0 60])
set(gca, 'xtick', [1 : length(wristExcursion)], 'xticklabel', freqs)
title('Elbow excursion')

subplot(1, 3, 2)
hold on
for iF = 1 : length(wristSlope)
    plot(iF*ones(size(wristSlope{iF}))+ 0.2*(rand(length(wristSlope{iF}),1)-0.5), wristSlope{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(wristSlope{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

%ylim([0 2])
set(gca, 'xtick', [1 : length(wristExcursion)], 'xticklabel', freqs)
title('Elbow excursion slope')

subplot(1, 3, 3)
hold on
for iF = 1 : length(wristMaxDer)
    plot(iF*ones(size(wristMaxDer{iF}))+ 0.2*(rand(length(wristSlope{iF}),1)-0.5), wristMaxDer{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(wristMaxDer{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end
ylim([0 12])
set(gca, 'xtick', [1 : length(wristExcursion)], 'xticklabel', freqs)
title('Peak Elbow velocity')



%% Figure on Wrist excursion

figure

subplot(1, 3, 1)
hold on
for iF = 1 : length(wristExcursion)
    find(isoutlier(wristExcursion{iF}))
    plot(iF*ones(size(wristExcursion{iF}))+ 0.2*(rand(length(wristSlope{iF}),1)-0.5), wristExcursion{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(wristExcursion{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

ylim([0 35])
set(gca, 'xtick', [1 : length(wristExcursion)], 'xticklabel', freqs)
title('Wrist excursion')

subplot(1, 3, 2)
hold on
for iF = 1 : length(wristSlope)
    plot(iF*ones(size(wristSlope{iF}))+ 0.2*(rand(length(wristSlope{iF}),1)-0.5), wristSlope{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(wristSlope{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

%ylim([0 2])
set(gca, 'xtick', [1 : length(wristExcursion)], 'xticklabel', freqs)
title('Wrist excursion slope')

subplot(1, 3, 3)
hold on
for iF = 1 : length(wristMaxDer)
    plot(iF*ones(size(wristMaxDer{iF}))+ 0.2*(rand(length(wristSlope{iF}),1)-0.5), wristMaxDer{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(wristMaxDer{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end
ylim([0 12])
set(gca, 'xtick', [1 : length(wristExcursion)], 'xticklabel', freqs)
title('Peak Wrist velocity')




%% Figure on pinkie excursion

figure

subplot(1, 3, 1)
hold on
for iF = 1 : length(pinkieExcursion)
    plot(iF*ones(size(pinkieExcursion{iF}))+ 0.2*(rand(length(wristSlope{iF}),1)-0.5), abs(pinkieExcursion{iF}), 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(-pinkieExcursion{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

%ylim([0 40])
set(gca, 'xtick', [1 : length(pinkieExcursion)], 'xticklabel', freqs)
title('Pinkie excursion')

subplot(1, 3, 2)
hold on
for iF = 1 : length(pinkieSlope)
    plot(iF*ones(size(pinkieSlope{iF}))+ 0.2*(rand(length(pinkieSlope{iF}),1)-0.5), pinkieSlope{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(pinkieSlope{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

%ylim([0 2])
set(gca, 'xtick', [1 : length(pinkieExcursion)], 'xticklabel', freqs)
title('Pinkie excursion slope')

subplot(1, 3, 3)
hold on
for iF = 1 : length(pinkieMaxDer)
    plot(iF*ones(size(pinkieMaxDer{iF}))+ 0.2*(rand(length(pinkieSlope{iF}),1)-0.5), pinkieMaxDer{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(pinkieMaxDer{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end
%ylim([0 12])
set(gca, 'xtick', [1 : length(pinkieExcursion)], 'xticklabel', freqs)
title('Peak Pinkie velocity')
