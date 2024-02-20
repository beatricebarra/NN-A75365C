expDate = '20190702';
pinNames = {'E1', 'E2', 'E3', 'E4', 'E5', 'E6'};
muscles = {'DEL'; 'BIC'; 'TRI'; 'EDC'; 'FCR'; 'ECR'; 'FDS'; 'ABP'};
fS = 12207.03;
fsTDT = fS;
recsystems = {'TDT'};
if pittPC ==1
    Root = fullfile('R:\data_raw\primate\SCS_Monkey\', Animal, expDate, 'TERMINAL', 'BrienneTerminal');
else
    Root = fullfile('/Volumes/DATABACKUP1/SERVER/Data/Experiments', Animal, expDate, 'TERMINAL', 'BrienneTerminal');
end
filenameTDT = [expDate,'_BT', '-'];
filenameVICON = [expDate,'_',Animal, '_'];
filenameBKR = [expDate,'_',Animal, '_'];


selectedFiles = [26 : 31 ; ...%PIN5
                ]; 
selectedFilesList = selectedFiles';
selectedFilesList = selectedFilesList(:)';

pins = [5];
npins = size(selectedFiles, 1);
filesamp = [850*ones(1,6)];

amps = sort(unique(filesamp));        
filesfreqs = repmat([20 40 50 60 80 100], npins, 1)
freqs = [20 40 50 60 80 100 ];

