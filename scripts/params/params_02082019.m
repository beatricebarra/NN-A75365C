expDate = '20190802';

% Files from 5 to 7 have no kinematics
TDTfileNames = {'20190802_Ygritte-1', ...
    '20190802_Ygritte-2', ...
    '20190802_Ygritte-3',...
    '20190802_Ygritte-4',...
    };

VICONfileNames = {'20190802_Ygritte_01_M', ...
    '20190802_Ygritte_02_M', ...
    '20190802_Ygritte_03_M',...
    '20190802_Ygritte_04_M',...
    };

pinNames = {};
muscles = {'DEL'; 'BIC'; 'TRI'; 'EDC'; 'FCR'; 'ECR'; 'FDS'; 'ABP'};
currents= {
    };
fsTDT = 12207.03;
fsVICON_analog = 1000;
fsVICON_kin = 100;
fS = fsTDT;
if pittPC ==1
    %Root = fullfile('D:\SERVER\Data\Experiments', Animal, expDate);
    Root = fullfile('R:\data_raw\primate\SCS_Monkey\', Animal, expDate);
    currentpath = 'C:\GIT\unifr_cervicalEES\';
else
    Root = fullfile('/Volumes/DATABACKUP1/SERVER/Data/Experiments', Animal, expDate);
end

recsystems = {'TDT', 'VICON'};