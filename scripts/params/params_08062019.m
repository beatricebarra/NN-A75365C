expDate = '20190608';
TDTfileNames = {...
    '20190608_Brienne-5',...
    '20190608_Brienne-6'};

VICONfileNames = {'20190608_Brienne_05', ...
    '20190608_Brienne_06', ...
    };
BKRfileNames = { ...
    '20190608_Brienne_005', ...
    '20190608_Brienne_006',...
    };
pinNames = {};
selpins = [2 3 6];
muscles = {'DEL'; 'BIC'; 'TRI'; 'EDC'; 'FCR'; 'ECR'; 'FDS'; 'ABP'};
currents= {
    };
fsTDT = 12207.03;
fsVICON_analog = 1000;
fsVICON_kin = 100;
fS = fsTDT;
if pittPC ==1
    Root = fullfile('E:\SERVER\Data\Experiments', Animal, expDate);
else
    Root = fullfile('/Volumes/DATABACKUP1/SERVER/Data/Experiments', Animal, expDate);
end

recsystems = {'TDT', 'VICON', 'BKR'};
