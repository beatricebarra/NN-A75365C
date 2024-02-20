Animal = 'Brienne';
expDate = '20190619';
pinNames = {'E1', 'E2', 'E3', 'E4', 'E5', 'E6'};
muscles = {'DEL'; 'BIC'; 'TRI'; 'EDC'; 'FCR'; 'ECR'; 'FDS'; 'ABP'};
fS = 12207.03;
fsTDT = fS;
recsystems = {'TDT'};
if pittPC ==1
   
    Root = fullfile('R:\data_raw\primate\SCS_Monkey', Animal, expDate, 'TDT');
else
    Root = fullfile('/Volumes/DATABACKUP1/SERVER/Data/Experiments', Animal, expDate, 'TDT');
end
filenameTDT = [expDate,'_',Animal,'-'];
filenameVICON = [expDate,'_',Animal, '_'];
filenameBKR = [expDate,'_',Animal, '_'];


selectedFiles = [...
                 15 14 17 13 12 16; ...%PIN2
                 25 26 27 28 29 30; ...% PIN5
                ]; 
selectedFilesList = selectedFiles';
selectedFilesList = selectedFilesList(:)';

pins = [2 5];
npins = size(selectedFiles, 1);
filesamp = [360*ones(1,6); ...
            380*ones(1,6)];

amps = sort(unique(filesamp));        
filesfreqs = repmat([20 40 50 60 80 100], npins, 1)
freqs = [20 40 50 60 80 100 ];







