expDate = '20190814';
fileNames = {'20190814_Ygritte-1', ...
    '20190814_Ygritte-2', ...
    '20190814_Ygritte-3',...
    '20190814_Ygritte-4',...
    '20190814_Ygritte-5',...
    '20190814_Ygritte-6'};

pinNames = {'E1', 'E2', 'E3', 'E4', 'E5', 'E6'};
muscles = {'DEL'; 'BIC'; 'TRI'; 'EDC'; 'FCR'; 'ECR'; 'FDS'; 'ABP'};
currents= {[500:70:1200]; ...
    [100:90:1000]; ...
    [300:90:1200]; ...
    [100:80:900]; ...
    [700:130:2000]; ...
    [400:40:800]; ...
    };

if pittPC ==1
    Root = fullfile('R:\data_raw\primate\SCS_Monkey', Animal, expDate, 'TDT');
else
    Root = fullfile('/Volumes/DATABACKUP1/SERVER/Data/Experiments', Animal, expDate, 'TDT');
end