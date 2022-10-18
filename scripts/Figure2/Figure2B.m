
%% Rostral Index
datapath = 'pathtodata'; 

MkBr=  readtable(fullfile(datapath, 'Figure2B_E2_Mk-Br.csv')); 
MkYg=  readtable(fullfile(datapath, 'Figure2B_E4_Mk-Yg.csv')); 
MkSa=  readtable(fullfile(datapath, 'Figure2B_E3_Mk-Sa.csv')); 

dataset_rostral = {MkBr, MkYg, MkSa}
rostrocaudalindex = [8 7 6 3 4 5 2 ];
[sortrcindex, isort] = sort(rostrocaudalindex); 
nanimalsmuscles = [2 3 3 2 2 2 3 2];
colorsgrey = [125 125 125]./255;
ordered_muscle_set = {'DEL', 'BIC', 'TRI', 'FCR', 'ECR', 'EDC', 'FDS'};


labels = MkBr.Properties.VariableNames'; 
muscles = labels(2:end, 1); 
nChan = size(MkBr, 2) -2;
subplot(1, 2, 1)
hold on
meanAnimals = zeros(1, 8);
for iA = 1 : length(dataset_rostral)
    animal = dataset_rostral{iA}; 
    for ch = 1 : nChan
        x = animal.Current; 
        try
            y = animal.(labels{ch+1}); 
        catch
            y = []; 
        end
    
        activationindex(ch) = mean(y);
        plot(activationindex(ch), rostrocaudalindex(ch),  'o', 'Linewidth', 2,...
            'MarkerFaceColor', colorsgrey, ...
            'MarkerEdgeColor', colorsgrey, ...
            'Markersize', 10)
        if isnan(activationindex(ch))
            meanAnimals(1, ch) = meanAnimals(1, ch);
    
        else
            meanAnimals(1, ch) = meanAnimals(1, ch) + activationindex(ch);
        end
    
    end
end


meanAnimals = meanAnimals./nanimalsmuscles;
plot(meanAnimals(isort), sortrcindex, 'k')
plot(meanAnimals(isort), sortrcindex,  'o', 'Linewidth', 2,...
                'MarkerFaceColor', [0 0 0], ...
                'MarkerEdgeColor', [0 0 0], ...
                'Markersize', 16)
set(gca, 'yticklabel', fliplr(ordered_muscle_set))
title('C6/C7 ')

%% Caudal Index
datapath = 'pathtodata'; 

MkBr=  readtable(fullfile(datapath, 'Figure2B_E5_Mk-Br.csv')); 
MkYg=  readtable(fullfile(datapath, 'Figure2B_E6_Mk-Yg.csv')); 
MkSa=  readtable(fullfile(datapath, 'Figure2B_E4_Mk-Sa.csv')); 

dataset_rostral = {MkBr, MkYg, MkSa}
rostrocaudalindex = [8 7 6 3 4 5 2 ];
[sortrcindex, isort] = sort(rostrocaudalindex); 
nanimalsmuscles = [2 3 3 2 2 2 3 2];
colorsgrey = [125 125 125]./255;
ordered_muscle_set = {'DEL', 'BIC', 'TRI', 'FCR', 'ECR', 'EDC', 'FDS'};


labels = MkBr.Properties.VariableNames'; 
muscles = labels(2:end, 1); 
nChan = size(MkBr, 2) -2;
subplot(1, 2, 2)
hold on
meanAnimals = zeros(1, 8);
for iA = 1 : length(dataset_rostral)
    animal = dataset_rostral{iA}; 
    for ch = 1 : nChan
        x = animal.Current; 
        try
            y = animal.(labels{ch+1}); 
        catch
            y = []; 
        end
    
        activationindex(ch) = mean(y);
        plot(activationindex(ch), rostrocaudalindex(ch),  'o', 'Linewidth', 2,...
            'MarkerFaceColor', colorsgrey, ...
            'MarkerEdgeColor', colorsgrey, ...
            'Markersize', 10)
        if isnan(activationindex(ch))
            meanAnimals(1, ch) = meanAnimals(1, ch);
    
        else
            meanAnimals(1, ch) = meanAnimals(1, ch) + activationindex(ch);
        end
    
    end
end


meanAnimals = meanAnimals./nanimalsmuscles;
plot(meanAnimals(isort), sortrcindex, 'k')
plot(meanAnimals(isort), sortrcindex,  'o', 'Linewidth', 2,...
                'MarkerFaceColor', [0 0 0], ...
                'MarkerEdgeColor', [0 0 0], ...
                'Markersize', 16)
set(gca, 'yticklabel', fliplr(ordered_muscle_set))
title('T1')