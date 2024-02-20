Animal = 'Ygritte';
expDate = '20190814';
pinNames = {'E1', 'E2', 'E3', 'E4', 'E5', 'E6'};
muscles = {'DEL'; 'BIC'; 'TRI'; 'EDC'; 'FCR'; 'ECR'; 'FDS'; 'ABP'};
fS = 12207.03;
fsTDT = fS;
recsystems = {'TDT'};
if pittPC ==1
    Root = fullfile('R:\data_raw\primate\SCS_Monkey', Animal, expDate);
else
    Root = fullfile('/Volumes/DATABACKUP1/SERVER/Data/Experiments', Animal, expDate);
end
filenameTDT = [expDate,'_', Animal,'-'];
filenameVICON = [expDate,'_',Animal, '_'];
filenameBKR = [expDate,'_',Animal, '_'];


selectedFiles = [...
                 8 13 9 7 10 11 12; ... PIN1
                 14 15 16 17 19 20 21 ; ...%PIN2
                 40 41 42 43 44 45 46; ...% PIN3
                 26 27 28 29 30 31 32; ...% PIN4
                 33 34 35 36 37 38 39; ...% PIN6
                ]; 
selectedFilesList = selectedFiles';
selectedFilesList = selectedFilesList(:)';

pins = [1 2 3 4 6];
npins = size(selectedFiles, 1);
filesamp = [800*ones(1,7); ...
            700*ones(1,7); ...
            900*ones(1,7);...
            650*ones(1,7); ...
            600*ones(1,7)];
amps = sort(unique(filesamp));        
filesfreqs = repmat([20 30 40 50 80 100 120], npins, 1)
freqs = [20 30 40 50 80 100 120];

