clear all
close all
Sansa_fig = '/Users/Bea/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/REVIEW-ROUND2/Figure5/Mk_Sa_EMG_allmuscles.fig';

fig = openfig(Sansa_fig);
fig = gcf;
axesObjs = fig.Children



%% Create table


baseline_plot_idx = [ 1.5]; 
stim_plot_idx = [2.5]; 


muscle_idx = [1 2 3];
muscle_axes_idx = [7 6 2]; % BIC TRI FDS
muscle_labels = {'BIC', 'TRI', 'FDS'}; 
iMM = 0;
clear stim_std stim_mean bas_mean bas_std
for iM = muscle_axes_idx
    myax = axesObjs(iM);
    iMM = iMM + 1;
    elements = myax.Children; 
    Bars = findobj(elements, 'Type','Bar'); 
    ErrorBars = findobj(elements, 'Type','ErrorBar')
    
    for iB = 1 : length(Bars)
        if mod(floor(Bars(iB).XData), 2) % Odd number = baseline condition
            date_idx = find(baseline_plot_idx == Bars(iB).XData); 
            bas_mean{iMM}(date_idx) = Bars(iB).YData; 
        else % Even number = stim condition
            date_idx = find(stim_plot_idx == Bars(iB).XData); 
            stim_mean{iMM}(date_idx) = Bars(iB).YData; 
        end
        
        if mod(floor(ErrorBars(iB).XData), 2) % Odd number = baseline condition
            date_idx = find(baseline_plot_idx == ErrorBars(iB).XData); 
            bas_std{iMM}(date_idx) = ErrorBars(iB).YPositiveDelta; 
        else % Even number = stim condition
            date_idx = find(stim_plot_idx == ErrorBars(iB).XData); 
            stim_std{iMM}(date_idx) = ErrorBars(iB).YPositiveDelta; 
        end
    end
    
end



%% Write table 

mytable = cell(13,7);
mytable(:, 1) = {'';'BIC'; '';'';''; 'TRI'; '';'';''; 'FDS'; '';'';''}; 
mytable(:, 2) = {'';'Baseline'; ''; 'Stim'; ''; 'Baseline'; ''; 'Stim'; '';  'Baseline'; ''; 'Stim'; ''; }; 
mytable(:, 3) = {'';'Mean'; 'Std'; 'Mean'; 'Std'; 'Mean'; 'Std'; 'Mean'; 'Std'; 'Mean'; 'Std'; 'Mean'; 'Std'; }; 

for iM = 1 : 3
    mytable((iM-1)*4 + 2, 4:6) = num2cell(bas_mean{iM}./(10^-12)); 
    mytable((iM-1)*4 + 3, 4:6) = num2cell(bas_std{iM}./(10^-12)); 
    mytable((iM-1)*4 + 4, 4:6) = num2cell(stim_mean{iM}./(10^-12)); 
    mytable((iM-1)*4 + 5, 4:6) = num2cell(stim_std{iM}./(10^-12)); 
end

filename = '/Users/Bea/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/REVIEW-ROUND2/RawDataSheets/Figure5A.xlsx';
writecell(mytable,filename,'Sheet','Sansa')
