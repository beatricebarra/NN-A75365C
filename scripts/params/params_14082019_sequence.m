expDate = '20190814';
TDTfileNames = {'20190814_Ygritte-47', ...
    '20190814_Ygritte-48'};

pinNames = {'E1', 'E2', 'E3', 'E4', 'E5', 'E6'};
muscles = {'DEL'; 'BIC'; 'TRI'; 'EDC'; 'FCR'; 'ECR'; 'FDS'; 'ABP'};
fS = 12207.03;
fsTDT = fS;
recsystems = {'TDT'};
if pittPC ==1
    Root = fullfile('D:\SERVER\Data\Experiments', Animal, expDate);
else
    Root = fullfile('/Volumes/DATABACKUP1/SERVER/Data/Experiments', Animal, expDate);
end