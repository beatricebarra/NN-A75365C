expDate = '20190603';
fileNames = {'20190603_Brienne_-1', ...
    '20190603_Brienne_-3', ...
    '20190603_Brienne_-5',...
    '20190603_Brienne_-6',...
    '20190603_Brienne_-7',...
    '20190603_Brienne_-9'};

pinNames = {'E1', 'E2', 'E3', 'E4', 'E5', 'E6'};
muscles = {'DEL'; 'BIC'; 'TRI'; 'EDC'; 'FCR'; 'ECR'; 'FDS'; 'ABP'};
currents= {[100:80:720]; ...
    [200:50:700]; ...
    [100:50:600]; ...
    [100:50:600]; ...
    [100:60:700]; ...
    [500:50:1000]; ...
    };
if pittPC ==1
    Root = fullfile('P:\data_raw\primate\SCS_Monkey', Animal, expDate, 'TDT');
else
    Root = fullfile('/Volumes/DATABACKUP1/SERVER/Data/Experiments', Animal, expDate, 'TDT');
end