sansa_fig = '/Users/Bea/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/REVIEW-ROUND2/Figure4/final/MkSa_rates.fig';

fig = openfig(sansa_fig);
fig = gcf;
axesObjs = fig.Children

reachAx = axesObjs(3)
graspAx = axesObjs(2)
pullAx = axesObjs(1)


%% Build table 

mytable = cell(13,7);
mytable(:, 1) = {'';'Reach'; '';'';''; 'Grasp'; '';'';''; 'Pull'; '';'';''}; 
mytable(:, 2) = {'';'Baseline'; ''; 'Stim'; ''; 'Baseline'; ''; 'Stim'; '';  'Baseline'; ''; 'Stim'; ''; }; 
mytable(:, 3) = {'';'Mean'; 'Std'; 'Mean'; 'Std'; 'Mean'; 'Std'; 'Mean'; 'Std'; 'Mean'; 'Std'; 'Mean'; 'Std'; }; 


%% Reach

dataObjs = reachAx.Children
for i = length(dataObjs) : -1 : 2
    
    if dataObjs(i).Color(1) == 0
        if length(dataObjs(i).XData)<2
            mean_bas = dataObjs(i).YData;
        else
            std_bas(dataObjs(i).XData(1)) = (dataObjs(i).YData(2) - dataObjs(i).YData(1))/2; 
        end
    else
        if length(dataObjs(i).XData)<2
            mean_stim = dataObjs(i).YData;
        else
            std_stim(dataObjs(i).XData(1)) = (dataObjs(i).YData(2) - dataObjs(i).YData(1))/2; 
        end
    end
end


mytable(2, 4:7) = num2cell(60*mean_bas); 
mytable(3, 4:7) = num2cell(60*std_bas); 
mytable(4, 4:7) = num2cell(60*mean_stim); 
mytable(5, 4:7) = num2cell(60*std_stim); 



filename = '/Users/Bea/Dropbox/CERVICAL SCS PAPER/NATURE NEURO/REVIEW-ROUND2/RawDataSheets/Figure4C.xlsx';
writecell(mytable,filename,'Sheet','Sansa')
