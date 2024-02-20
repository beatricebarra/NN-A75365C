%DLCanalysis

clear all
DLCfiles_folder = 'P:\data_raw\primate\SCS_Monkey\Brienne\20190702\TERMINAL\labeled_videos';
tracking_folder = 'P:\data_raw\primate\SCS_Monkey\Brienne\20190702\TERMINAL\trackedfiles';
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


%FLexion set 
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


figure
cols = ['r', 'b', 'g', 'y', 'k', 'm'];
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
    if strcmp(seltrackedfiles{1}, 'C0294_reachEvents.mat')
        if iF == 6
             kin.starts(8) = [];
            kin.starts(end) = [];
        end
    end
        
    cut = min(kin.ends-kin.starts);
    
    filtered_wrist = lowPass(kin.angles.wrist_flex', thisfs, 6, 3);
    filtered_pinkie = lowPass(kin.angles.pinkie_flex', thisfs, 6, 3);
    
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
        halfstimtime = floor(length(wristsnip)/2)
        [maxpeak, maxpeakloc]= max(wristsnip(1 : floor(movlength/2)));
        %plot(wristsnip, 'Color', cols(iF))
        %hold on
        %plot(wristsnip)
        %hold on, plot(maxpeakloc, maxpeak, 'o')
        
        
        wristSnips{iF}(:,imov) = wristsnip(1:cut) - wristsnip(1);
        wristSlope{iF}(imov) = (maxpeak - filtered_wrist(kin.starts(imov)))/(maxpeakloc);
        wristExcursion{iF}(imov) = wristsnip(halfstimtime)- wristsnip(1);
        wristRiseTime{iF}(imov) =(maxpeakloc);
        wristMaxDer{iF}(imov) =max(diff(wristsnip));
        % Pinkie
        pinkiesnip = filtered_pinkie(kin.starts(imov): kin.ends(imov));
        plot(pinkiesnip, 'Color', cols(iF))
        hold on
        
        [minpeak, minpeakloc]= min(pinkiesnip(1 : floor(movlength/2)));

        pinkieSnips{iF}(:,imov) = pinkiesnip(1:cut) - pinkiesnip(1);
        pinkieSlope{iF}(imov) = (minpeak - filtered_wrist(kin.starts(imov)))/(minpeakloc);
        %pinkieExcursion{iF}(imov) = pinkiesnip(minpeakloc)- pinkiesnip(1);%
        pinkieExcursion{iF}(imov) = pinkiesnip(1) - pinkiesnip(minpeakloc) ;% Ygritte version
        pinkieRiseTime{iF}(imov) =(minpeakloc);
        pinkieMaxDer{iF}(imov) =max(diff(pinkiesnip));
        
        
        
    end
end




%% Figure on Elbow. wrist and pinkie excursion

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

% Figure on Wrist excursion
subplot(1, 3, 2)
hold on
for iF = 1 : length(wristExcursion)
    find(isoutlier(wristExcursion{iF}))
    plot(iF*ones(size(wristExcursion{iF}))+ 0.2*(rand(length(wristExcursion{iF}),1)-0.5), wristExcursion{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(wristExcursion{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end


set(gca, 'xtick', [1 : length(wristExcursion)], 'xticklabel', freqs)
title('Wrist excursion')


% Figure on Pinkie excursion
subplot(1, 3, 3)
hold on
for iF = 1 : length(pinkieExcursion)
    find(isoutlier(pinkieExcursion{iF}))
    plot(iF*ones(size(pinkieExcursion{iF}))+ 0.2*(rand(length(pinkieExcursion{iF}),1)-0.5), pinkieExcursion{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(pinkieExcursion{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

%ylim([0 35])
set(gca, 'xtick', [1 : length(pinkieExcursion)], 'xticklabel', freqs)
title('Pinkie excursion')


%% Statistics Elbow
data_vect = [];
group_vect = [];
for iF = 1 : length(elbowExcursion)
    data_vect = [data_vect, elbowExcursion{iF}];
    group_vect = [group_vect, iF*ones(size(elbowExcursion{iF}))]; 
    data_cell{iF} = elbowExcursion{iF};
end
ncomparisons = 5;
ref1 = 0.05/ncomparisons;
ref2 = 0.01/ncomparisons; 
ref3 = 0.001/ncomparisons;

disp('Group 1 vs Group 2')
p = ranksum(data_cell{1}, data_cell{2}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 1 vs Group 3')
p = ranksum(data_cell{1}, data_cell{3}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 1 vs Group 4')
p = ranksum(data_cell{1}, data_cell{4}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 1 vs Group 5')
p = ranksum(data_cell{1}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 1 vs Group 6')
p = ranksum(data_cell{1}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

%[P,ANOVATAB,STATS] = kruskalwallis(data_vect, group_vect); 
%multcompare(STATS)
%% Statistics Wrist
data_vect = [];
group_vect = [];
for iF = 1 : length(wristExcursion)
    data_vect = [data_vect, wristExcursion{iF}];
    group_vect = [group_vect, iF*ones(size(wristExcursion{iF}))];
    data_cell{iF} = wristExcursion{iF};
end
%[P,ANOVATAB,STATS] = kruskalwallis(data_vect, group_vect); 
%multcompare(STATS)
ncomparisons = 5;
ref1 = 0.05/ncomparisons;
ref2 = 0.01/ncomparisons; 
ref3 = 0.001/ncomparisons;

disp('Group 1 vs Group 2')
p = ranksum(data_cell{1}, data_cell{2}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 2 vs Group 3')
p = ranksum(data_cell{2}, data_cell{3}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 3 vs Group 4')
p = ranksum(data_cell{3}, data_cell{4}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 4 vs Group 5')
p = ranksum(data_cell{4}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 5 vs Group 6')
p = ranksum(data_cell{5}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])


%% Statistics Fingers
data_vect = [];
group_vect = [];
for iF = 1 : length(pinkieExcursion)
    data_vect = [data_vect, pinkieExcursion{iF}];
    group_vect = [group_vect, iF*ones(size(pinkieExcursion{iF}))]; 
end

[P,ANOVATAB,STATS] = kruskalwallis(data_vect, group_vect); 
multcompare(STATS)

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
