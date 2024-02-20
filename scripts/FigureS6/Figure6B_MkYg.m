%DLCanalysis
clear all
DLCfiles_folder = 'P:\data_generated\primate\SCS_Monkey\Ygritte\20190814\DLC\kinematics\';
tracking_folder = 'P:\data_generated\primate\SCS_Monkey\Ygritte\20190814\DLC\trackedfiles\';


% Flex set 
selDLCfiles = {'C0476RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0477RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0478RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0479RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0480RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0481RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0482RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    };


seltrackedfiles = {'C0476_reachEvents.mat'; ...
    'C0477_reachEvents.mat'; ...
    'C0478_reachEvents.mat'; ...
    'C0479_reachEvents.mat'; ...
    'C0480_reachEvents.mat'; ...
    'C0481_reachEvents.mat'; ...
    'C0482_reachEvents.mat'; ...
    };


% %Ext set 
selDLCfiles = {'C0483RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0484RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0485RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0486RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0487RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0488RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    'C0489RCDLC_resnet50_An_MkYgMay11shuffle1_1030000.csv'; ...
    };


seltrackedfiles = {'C0483_reachEvents.mat'; ...
    'C0484_reachEvents.mat'; ...
    'C0485_reachEvents.mat'; ...
    'C0486_reachEvents.mat'; ...
    'C0487_reachEvents.mat'; ...
    'C0488_reachEvents.mat'; ...
    'C0489_reachEvents.mat'; ...
    };

fs = 60;
freqs = [20 30 40 50 80 100 120  ];

freqsLab = {'20' '30' '40' '50' '80' '100' '120'};


figure
hold on
cols = ['r', 'b', 'g', 'y', 'k', 'm', 'r'];
for iF = 1 : length(selDLCfiles)
    
    selDLCfile = fullfile(DLCfiles_folder, selDLCfiles{iF});
    seltrackedfile =  fullfile(tracking_folder, seltrackedfiles{iF});
    kin = importDLCkinematics(selDLCfile, 'Ygritte', '20190814')
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
        elbowVelocity{iF}(:,imov) = diff(elbowSnips{iF}(:,imov) );
        elbowRiseTime{iF}(imov) =(minpeakloc);
        
        elbowMaxDer{iF}(imov) =max(diff(elbowsnip));
        
        % Wrist
        wristsnip = filtered_wrist(kin.starts(imov): kin.ends(imov));
        [maxpeak, maxpeakloc]= max(wristsnip);
        halfstimtime = floor(length(wristsnip)/2)
         plot(wristsnip, 'Color', cols(iF))
         ylim([0 60])
        hold on
        wristSnips{iF}(:,imov) = wristsnip(1:cut) - wristsnip(1);
        wristSlope{iF}(imov) = (maxpeak - filtered_wrist(kin.starts(imov)))/(maxpeakloc);
        wristExcursion{iF}(imov) = wristsnip(halfstimtime)- wristsnip(1);
        wristVelocity{iF}(:,imov) = diff(wristSnips{iF}(:, imov) );
        wristRiseTime{iF}(imov) =(maxpeakloc);
        wristMaxDer{iF}(imov) =max(diff(wristsnip));
        
        % Pinkie
        pinkiesnip = filtered_pinkie(kin.starts(imov): kin.ends(imov));
        [minpeak, minpeakloc]= min(pinkiesnip(1 : floor(movlength/2)));
        pinkieSnips{iF}(:,imov) = pinkiesnip(1:cut) - pinkiesnip(1);
        pinkieSlope{iF}(imov) = (minpeak - filtered_wrist(kin.starts(imov)))/(minpeakloc);
        pinkieExcursion{iF}(imov) = pinkiesnip(1) - pinkiesnip(minpeakloc) ;
        pinkieVelocity{iF}(:, imov) = diff(pinkieSnips{iF}(:, imov) );
        pinkieRiseTime{iF}(imov) =(minpeakloc);
        pinkieMaxDer{iF}(imov) =max(diff(pinkiesnip));
        %plot(pinkiesnip, 'Color', cols(iF))
         %ylim([0 60])
        %hold on
        
        % Smoothness
       
        xwrist = kin.wrist_l.x(kin.starts(imov): kin.ends(imov));
        zwrist = kin.wrist_l.y(kin.starts(imov): kin.ends(imov));
        lengthMOV = sum(sqrt(diff(xwrist).^2 + diff(zwrist).^2 ));
        jerkx = diff(diff(diff(xwrist))); jerkz = diff(diff(diff(zwrist))); 
        jerk = sqrt(jerkx.^2+ jerkz.^2);
        durMOV = (kin.ends(imov)- kin.starts(imov))*1/fs;
        smoothness{iF}(imov) = sqrt(1/3*sum(abs(jerk))*(durMOV.^5/lengthMOV.^2));
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

ylim([-10 60])
set(gca, 'xtick', [1 : length(wristExcursion)], 'xticklabel', freqs)
title('Wrist excursion')


% Figure on Pinkie excursion
subplot(1, 3, 3)
hold on
for iF = 1 : length(pinkieExcursion)
    iF
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
ncomparisons = 8;
ref1 = 0.05/ncomparisons;
ref2 = 0.01/ncomparisons; 
ref3 = 0.001/ncomparisons;

disp('Group 1 vs Group 5')
p = ranksum(data_cell{1}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 1 vs Group 6')
p = ranksum(data_cell{1}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 2 vs Group 5')
p = ranksum(data_cell{2}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 2 vs Group 6')
p = ranksum(data_cell{2}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 3 vs Group 5')
p = ranksum(data_cell{3}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 3 vs Group 6')
p = ranksum(data_cell{3}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 4 vs Group 5')
p = ranksum(data_cell{4}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 4 vs Group 6')
p = ranksum(data_cell{4}, data_cell{6}); 
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
ncomparisons = 8;
ref1 = 0.05/ncomparisons;
ref2 = 0.01/ncomparisons; 
ref3 = 0.001/ncomparisons;

disp('Group 1 vs Group 5')
p = ranksum(data_cell{1}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 1 vs Group 6')
p = ranksum(data_cell{1}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 2 vs Group 5')
p = ranksum(data_cell{2}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 2 vs Group 6')
p = ranksum(data_cell{2}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 3 vs Group 5')
p = ranksum(data_cell{3}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 3 vs Group 6')
p = ranksum(data_cell{3}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 4 vs Group 5')
p = ranksum(data_cell{4}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 4 vs Group 6')
p = ranksum(data_cell{4}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

%% Statistics Fingers
data_vect = [];
group_vect = [];
for iF = 1 : length(pinkieExcursion)
    data_vect = [data_vect, pinkieExcursion{iF}];
    group_vect = [group_vect, iF*ones(size(pinkieExcursion{iF}))]; 
    data_cell{iF} =  pinkieExcursion{iF};
end

%[P,ANOVATAB,STATS] = kruskalwallis(data_vect, group_vect); 
%multcompare(STATS)
ncomparisons = 6;
ref1 = 0.05/ncomparisons;
ref2 = 0.01/ncomparisons; 
ref3 = 0.001/ncomparisons;

% disp('Group 1 vs Group 2')
% p = ranksum(data_cell{1}, data_cell{2}); 
% s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
% disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])
% 
% disp('Group 2 vs Group 3')
% p = ranksum(data_cell{2}, data_cell{3}); 
% s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
% disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])
% 
% disp('Group 3 vs Group 4')
% p = ranksum(data_cell{3}, data_cell{4}); 
% s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
% disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])
% 
% disp('Group 4 vs Group 5')
% p = ranksum(data_cell{4}, data_cell{5}); 
% s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
% disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])
% 
% disp('Group 5 vs Group 6')
% p = ranksum(data_cell{5}, data_cell{6}); 
% s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
% disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 1 vs Group 5')
p = ranksum(data_cell{1}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 1 vs Group 6')
p = ranksum(data_cell{1}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 2 vs Group 5')
p = ranksum(data_cell{2}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 2 vs Group 6')
p = ranksum(data_cell{2}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 3 vs Group 5')
p = ranksum(data_cell{3}, data_cell{5}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

disp('Group 3 vs Group 6')
p = ranksum(data_cell{3}, data_cell{6}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])


disp('Group 3 vs Group 4')
p = ranksum(data_cell{3}, data_cell{4}); 
s1 = p < ref1; s2 = p < ref2; s3 = p < ref3; 
disp(['pvalue = ', num2str(p), ': [' num2str(s1) num2str(s2) num2str(s3), ']'])

%% Figure on Elbow excursion (for extension dataset)

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



%% Figure on Wrist excursion (for flexion dataset)

figure

subplot(1, 3, 1)
hold on
for iF = 1 : length(wristExcursion)
    plot(iF*ones(size(wristExcursion{iF}))+ 0.2*(rand(length(wristSlope{iF}),1)-0.5), wristExcursion{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(wristExcursion{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

%ylim([0 40])
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




%% Figure on pinkie excursion (for flexion dataset)

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


%% Smoothness (Figure 4A)


figure
hold on
for iF = 1 : length(smoothness)
    plot(iF*ones(size(smoothness{iF}))+ 0.2*(rand(length(smoothness{iF}),1)-0.5), smoothness{iF}, 'o',...
        'MarkerFaceColor', [0.5, 0.5, 0.5], ...
        'MarkerEdgeColor', [0, 0, 0])
    plot(iF, median(smoothness{iF}), 'o',...
        'MarkerFaceColor', [0, 0, 0], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'MarkerSize', 12)
    
end

set(gca, 'xtick', [1 : length(smoothness)], 'xticklabel', freqs)
title('Peak Pinkie velocity')

