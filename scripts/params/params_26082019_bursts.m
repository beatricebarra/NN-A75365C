expDate = '20190826';
pinNames = {'E4'};
muscles = {'DEL'; 'BIC'; 'TRI'; 'EDC'; 'FCR'; 'ECR'; 'FDS'; 'ABP'};
fS = 12207.03;
fsTDT = fS;
recsystems = {'TDT'};
if pittPC ==1
    Root = fullfile('D:\SERVER\Data\Experiments', Animal, expDate);
else
    Root = fullfile('/Volumes/DATABACKUP1/SERVER/Data/Experiments', Animal, expDate);
end
filenameTDT = [expDate,'_', Animal,'-'];
filenameVICON = [expDate,'_',Animal, '_'];
filenameBKR = [expDate,'_',Animal, '_'];


selectedFiles = [...
                 31 32 33 34;... % PIN4
                ]; 
selectedFilesList = selectedFiles';
selectedFilesList = selectedFilesList(:)';

pins = [4];
npins = size(selectedFiles, 1);
filesamp = [800*ones(1,7); ...
           ]; 
amps = sort(unique(filesamp));        
freqs = [20 50 80 100];
filesfreqs = repmat(freqs, npins, 1)


